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


def get_df_mask(df, mask): 
    """
    
    Parameters
    ----------
    df : DataFrame
        DF with intervals
    mask : Set[Tuple[int, str]]
        Mask to apply. {(pos, aa)}

    Returns
    -------
    f_mask : Series
        combined mask

    """
    f_mask = True
    for i in [df[k] == v.upper() for k, v in mask]:
        f_mask &= i
    return f_mask


def get_motif(fg_intervals, bg_intervals, mask, current_pval, pval_threshold, count_theshold):
    """

    Parameters
    ----------
    fg_intervals : DataFrame
        DF with intervals.
    bg_intervals : DataFrame
        DF with intervals.
    mask : Set[Tuple[int, str]]
        Mask for motif. {(pos, aa)}
    current_pval : float
        current p_value for motif
    pval_threshold : float
        p_value threshold
    count_theshold : int
        minimum occurence to be considered for motif.

    Returns
    -------
    ans : List[List]
        List with motifs.[[mask, fg_size, bg_size, current_pval],..]

    """
    fg_slice = fg_intervals[get_df_mask(fg_intervals, mask)]
    bg_slice = bg_intervals[get_df_mask(bg_intervals, mask)]
    fg_size, bg_size = fg_slice.shape[0], bg_slice.shape[0]
    logging.debug(f'{mask=}, {current_pval=}, {pval_threshold=}')#' \n {fg_slice} \n {bg_slice}')
    fg_occ, bg_occ = utils.get_occurences(fg_slice), utils.get_occurences(bg_slice)
    fg_occ = fg_occ[fg_occ >= count_theshold]
    bg_prob = bg_occ / bg_size
    binom_prob = pd.DataFrame(binom.sf(fg_occ - 1, fg_size, bg_prob),
                              columns=bg_prob.columns, index=bg_prob.index)
    logging.debug(f"Sizes {fg_size=}, {bg_size=}")
    min_aa, min_pos = binom_prob.stack().idxmin()
    logging.debug(f"Candidate {min_aa=}, {min_pos=}", binom_prob.at[min_aa, min_pos])
    next_pval = binom_prob.at[min_aa, min_pos] * current_pval
    if binom_prob.at[min_aa, min_pos] < pval_threshold:
        next_mask = mask | {(min_pos, min_aa)}
        logging.debug(f'{mask=}, {next_mask=}, {fg_size}')
        ans = get_motif(fg_slice, bg_slice, next_mask, next_pval, pval_threshold, count_theshold)
    else:
        ans = [mask, fg_size, bg_size, current_pval]  
    return ans


def form_motif(mask, length):
    """

    Parameters
    ----------
    mask : Set[Tuple[int, str]]
        Mask for motif. {(pos, aa)}
    length : int
        motif length = length * 2 + 1

    Returns
    -------
    str
        Motif.
        
    """
    out = ['.'] * (2 * length + 1)
    for (pos, aa) in mask:
        if pos == 0:
            aa = aa.lower()
        out[pos + length] = aa
    return ''.join(out)


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
    DataFrame 
        Df with all found motifs.

    """
    current_pval = 1
    pval_threshold = args.p_value
    count_theshold = args.occurrences
    motifs = []
    initial_mask = {(0, args.modification_site)}
    while True:
        motif = get_motif(fg_intervals, bg_intervals, initial_mask, current_pval, pval_threshold, count_theshold)
        motif.append(fg_intervals.shape[0])
        motif.append(bg_intervals.shape[0])
        logging.debug("Motif found! ", motif)
        logging.debug("=====================")
        if motif[0] != initial_mask:
            motifs.append(motif)
            fg_intervals = fg_intervals[~get_df_mask(fg_intervals, motif[0])]
            bg_intervals = bg_intervals[~get_df_mask(bg_intervals, motif[0])]
        else:
            break
        if fg_intervals.empty:
            break
    out = pd.DataFrame(motifs, columns=['motif', 'fg_matches', 'bg_matches', 'p_value', 'fg_size', 'bg_size'])
    out['motif'] = out['motif'].apply(form_motif, args=[args.interval_length])
    return out.loc[:,['motif', 'p_value', 'fg_matches', 'fg_size', 'bg_size', 'bg_matches']]

