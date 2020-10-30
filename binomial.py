    # -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:55:43 2020

@author: Лиля
"""
import logging
import utils

import numpy as np
import pandas as pd

from pyteomics import fasta
from scipy.stats import binom


def get_binom_prob(fg_occ, bg_occ):
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
    fg_size = fg_occ[0].sum()
    bg_prob = bg_occ / bg_occ[0].sum()#х
    logging.debug("Foreground occ\n%s", fg_occ)
    logging.debug("Background probabilities\n%s", bg_prob)
    binom_prob = pd.DataFrame(binom.sf(fg_occ - 1, fg_size, bg_prob), columns=bg_prob.columns, index=bg_prob.index)
    logging.debug("Binomial probabilities\n%s", binom_prob)
    return binom_prob


def find_motifs(binom_prob, fg_occ, mask, args):
    """
    Parameters
    ----------
    binom_prob : TYPE
        binomial probabilities of aa occurences.
    fg_occ : 2-d array
        matrix with aa occurences among foreground intervals
    mask : str
        motif mask
    args : ArgumentParser
        ArgumentParset with parameters

    Returns
    -------
    motifs : List[List]
        List with [motif, p_value, counts] structure
    motifs_pairs: List[Tuple]
        List of motif tuples.

    """
    aa_df = pd.DataFrame(np.array([np.array(fg_occ.index)] * fg_occ.shape[1]).T, index=fg_occ.index, columns=fg_occ.columns)
    motifs = []
    motifs_pairs = []
    mot_length = len(mask)
    df_slice = aa_df[(binom_prob <= args.p_value) & (fg_occ >= args.occurrences)].apply(lambda x: x.dropna().to_dict(), axis=1)
    logging.debug("AA motif table:\n%s", df_slice)
    for i in df_slice:
        if i:
            for pos, aa in i.items():
                if aa != '-' and pos != 0:
                    motif = mask[:]
                    motifs_pairs.append((pos, aa))
                    index = pos + args.interval_length
                    if motif[index] == '.':
                        motif[index] = aa
                        mot = (''.join(motif).strip('.'), binom_prob.at[aa, pos], fg_occ.at[aa, pos])
                        if mot not in motifs:
                            motifs.append(mot)
    logging.debug("Motifs %s, mask %s", pd.DataFrame(motifs), mask)
    logging.debug("Motifs pairs %s", motifs_pairs)
    return motifs, motifs_pairs


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
    for m_motif in multi_list:
        for (pos, aa) in single_motif_pairs:
            if pos not in [i[0] for i in m_motif]:
                motif = m_motif + ((pos, aa),)
                if motif not in motifs_next:
                    motifs_next.append(motif)
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
    args : dict
        dict with parameters.

    Returns
    -------
    List with motifs, Dict with motifs.

    """
    motif_l = args.interval_length * 2 + 1
    print(f'{mask=}')
    logging.debug('FG %s mask is used\n%s', mask, fg_intervals[get_df_mask(fg_intervals, mask)] )
    fg_occ = utils.get_occurences(fg_intervals[get_df_mask(fg_intervals, mask)])
    if fg_occ.sum().sum() > 0:
        logging.debug("Occurance fg matrix:\n%s", fg_occ)
        logging.debug('BG %s mask is used\n%s', mask, bg_intervals[get_df_mask(bg_intervals, mask)] )
        bg_occ = utils.get_occurences(bg_intervals[get_df_mask(bg_intervals, mask)])
        logging.debug("Occurance bg matrix:\n%s", bg_occ)
        binom_prob = get_binom_prob(fg_occ, bg_occ)
        motif = ['.'] * motif_l
        for (pos, aa) in mask:
            motif[pos + args.interval_length] = aa
        logging.debug('Prepared mask %s', motif)
        motifs, motif_pairs = find_motifs(binom_prob, fg_occ, motif, args)
        return motifs, motif_pairs
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
    DataFrame with foreground intervals..
    args : dict
        Parameters

    Returns
    -------
    DataFrame with all found motifs.

    """
    final_motifs = []
    step = 1
    motifs, motifs_pairs = process_binom_stage(fg_intervals, bg_intervals, ((0, args.modification_site.lower()), ), args)
    final_motifs.extend(motifs)
    multi_motifs = get_multimotifs([((0, args.modification_site.lower()), ), ], motifs_pairs)
    logging.info(f'{step = }, new motifs {pd.DataFrame(motifs)} ')
    step += 1
    logging.debug('Double motifs candidates %s', multi_motifs)
    while True:
        tmp_mods = set()
        new_step_motifs = []
        for motif in multi_motifs:
            new_motifs, new_motifs_pairs = process_binom_stage(fg_intervals, bg_intervals, motif, args)
            logging.debug(f'{new_motifs_pairs=} \n {new_motifs=}')
            new_step_motifs.extend(new_motifs)
            tmp_mods.update(new_motifs_pairs)
        logging.info(f'{step = }, new motifs {pd.DataFrame(new_step_motifs)} \n {tmp_mods=}')
        final_motifs.extend(new_step_motifs)
        step += 1
        if tmp_mods:
            multi_motifs = get_multimotifs(multi_motifs, tmp_mods)
        else:
            break
    return pd.DataFrame(final_motifs)
    
