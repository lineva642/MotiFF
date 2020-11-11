    # -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:55:43 2020

@author: Лиля
"""
import logging
import utils

import numpy as np
import pandas as pd

from scipy.stats import binom


def get_binom_prob(fg_occ, bg_occ, bg_size, fg_size):
    """
    Parameters
    ----------
    fg_occ : DataFrame
        matrix with aa occurences among foreground intervals
    bg_occ : DataFrame
        matrix with aa occurences among background peptides
    Returns
    -------
    binom_prob : 2-d array
        Calculated binomial probabilities of aa occurences.
    """
    bg_prob = bg_occ / bg_size
    # logging.debug("Foreground occ\n%s", fg_occ)
    # logging.debug("Background probabilities\n%s", bg_prob)
    binom_prob = pd.DataFrame(binom.sf(fg_occ - 1, fg_size, bg_prob), columns=bg_prob.columns, index=bg_prob.index)
    # logging.debug("Binomial probabilities\n%s", binom_prob)
    return binom_prob


def find_motifs(binom_prob, fg_occ, mask, args):
    """
    Parameters
    ----------
    binom_prob : TYPE
        binomial probabilities of aa occurences.
    fg_occ : 2-d array
        matrix with aa occurences among foreground intervals
    mask : list[ tuples]
        motif mask [(pos, aa), (pos, aa), ]
    args : ArgumentParser
        ArgumentParset with parameters

    Returns
    -------
    motifs : List[List]
        List with [motif, p_value, counts] structure
    new_masks: List[List[Tuple]]
        List of new masks. [[(pos, aa), (pos, aa)], [(pos, aa)], ]

    """
    aa_df = pd.DataFrame(np.array([np.array(fg_occ.index)] * fg_occ.shape[1]).T, index=fg_occ.index, columns=fg_occ.columns)
    motifs = []
    new_masks = []
    df_slice = aa_df[(binom_prob <= args.p_value) & (fg_occ >= args.occurrences)].apply(lambda x: x.dropna().to_dict(), axis=1)
    motif = ['.'] * (2 * args.interval_length + 1)
    for (pos, aa) in mask:
        motif[pos + args.interval_length] = aa
    print(motif)
    for i in df_slice:
        if i:
            for pos, aa in i.items():
                index = pos + args.interval_length
                new_motif = motif[:]
                # print(f'{motif=}')
                if aa != '-' and motif[index] == '.':
                
                    new_masks.append(mask + [(pos, aa)])
                    new_motif[index] = aa
                    mot = (''.join(new_motif).strip('.'), binom_prob.at[aa, pos], fg_occ.at[aa, pos])
                    motifs.append(mot)
                    # print(mot, new_masks)
    # logging.debug("Motifs %s, mask %s", pd.DataFrame(motifs), mask)
    # logging.debug("Motifs pairs %s", new_masks)
    print(new_masks)
    return motifs, new_masks


def get_df_mask(df, mask): 
    """
    Parameters
    ----------
    df : DataFrame
        DF with intervals
    mask : tuple
        masks to apply

    Returns
    -------
    f_mask : Series
        combined mask

    """
    f_mask = True
    logging.debug(f'{mask=}')
    for i in [df[k] == v.upper() for k, v in mask]:
        f_mask &= i
    return f_mask

 
def get_multimotifs(multi_list, single_motif_pairs):
    """
    Prepares motif candidates for next step.

    Parameters
    ----------
    multi_list : List[Tuple]
        List of Tuples with motifs from previous step. [((pos, aa), (pos, aa))]
    single_motif_pairs : List[Tuple]
        List of Tuples with single aa motifs. [(pos, aa), (pos, aa)]

    Returns
    -------
    motifs_next : List[Tuple] of next motif candidates.

    """
    motifs_next = []
    print(len(single_motif_pairs), len(multi_list))
    for m_motif in multi_list:
        for (pos, aa) in single_motif_pairs:
            if pos not in [i[0] for i in m_motif]:
                motif = m_motif + [(pos, aa),]
                if sorted(motif) not in motifs_next:
                    motifs_next.append(sorted(motif))
    return motifs_next


def process_binom_stage(fg_intervals, bg_intervals, mask, args):
    """
    Finds motifs using Binomial probabilities.
    Parameters
    ----------
    fg_intervals : DataFrame
    DF with foreground intervals.
    bg_intervals : DataFrame
    DataFrame with foreground intervals.
    mask : Tuple[Tuple[pos, aa]]
        Motif mask. Tuple with fixed aa.
    args : ArgumentParser
        Parameters.

    Returns
    -------
    List with motifs, Dict with motifs.

    """
    motif_l = args.interval_length * 2 + 1
    # logging.debug('FG %s mask is used\n%s', mask, fg_intervals[get_df_mask(fg_intervals, mask)] )
    fg_occ = utils.get_occurences(fg_intervals[get_df_mask(fg_intervals, mask)])
    if fg_occ.sum().sum() > 0:
        
        # logging.debug('BG %s mask is used\n%s', mask, bg_intervals[get_df_mask(bg_intervals, mask)] )
        bg_occ = utils.get_occurences(bg_intervals[get_df_mask(bg_intervals, mask)])
        binom_prob = get_binom_prob(fg_occ, bg_occ, len(bg_intervals), len(fg_intervals))
        # logging.debug('Prepared mask %s', motif)
        motifs, new_masks = find_motifs(binom_prob, fg_occ, mask, args)
        return motifs, new_masks
    else:
        return [], {}


def binomial_alg(fg_intervals, bg_intervals, args):
    """
    Finds all motifs using Binomial probabilities.
    Parameters
    ----------
    fg_intervals : DataFrame
    DF with foreground intervals.
    bg_intervals : DataFrame
    DataFrame with background intervals.
    args : ArgumentParser
        Parameters
    Returns
    -------
    DataFrame with all found motifs.

    """
    final_motifs = []
    step = 1
    motif = [(0, args.modification_site.lower()),]
    motifs, masks = process_binom_stage(fg_intervals, bg_intervals, motif, args)
    final_motifs.extend(motifs)
    logging.info(f'{step = }, new motifs {pd.DataFrame(motifs)} ')
    step += 1
    while True:
        tmp_mods = []
        new_step_motifs = set()
        print(f'{step=}', len(masks))
        for motif in masks:
            new_motifs, new_masks = process_binom_stage(fg_intervals, bg_intervals, motif, args)
            new_step_motifs.update(new_motifs)
            tmp_mods.extend(new_masks)
            print(motif, new_masks)
            
        logging.info(f'{step = }, new motifs {pd.DataFrame(new_step_motifs)} \n')
        final_motifs.extend(list(new_step_motifs))
        step += 1
        if not new_masks:
            break
        else:
            masks = tmp_mods
    return pd.DataFrame(final_motifs)
    
