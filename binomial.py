# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:55:43 2020

@author: Лиля
"""

import os
from urllib.request import urlretrieve
import pylab
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyteomics import fasta
import gzip
#import seaborn as sns
from scipy.stats import binom
import math
import random
import profile
from scipy.stats import chisquare

import argparse
import logging
import utils


def get_binom_prob(fg_occ, bg_occ, args, results_saving_dir):
    """
    Parameters
    ----------
    fg_occ : 2-d array
        matrix with aa occurences among foreground intervals
    bg_occ : 2-d array
        matrix with aa occurences among background peptides
    args : ArgumentParser
        ArgumentParset with parameters
    results_saving_dir : str
        saving directory

    Returns
    -------
    binom_prob : 2-d array
        Calculated binomial probabilities of aa occurences.
    """
    
    fg_size = fg_occ[0].sum()
    bg_prob = bg_occ / bg_occ[0].sum()#х
    logging.debug("Background probabilities\n%s", bg_prob)
    binom_prob = pd.DataFrame(binom.sf(fg_occ - 1, fg_size, bg_prob), columns=bg_prob.columns, index=bg_prob.index)
    logging.debug("Binomial probabilities\n%s", binom_prob)
    utils.saving_table(results_saving_dir, binom_prob, args.interval_length, 'P_binomial_matrix')            
    logging.info(u'Binomial probability for each amino acid in matrix was counted')
    logging.debug("P_binomial matrix:\n%s", binom_prob)
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
    motifs : TYPE
        DESCRIPTION.
    motifs_dict : TYPE
        DESCRIPTION.

    """
    aa_df = pd.DataFrame(np.array([np.array(fg_occ.index)] * fg_occ.shape[1]).T, index=fg_occ.index, columns=fg_occ.columns)
    motifs = []
    motifs_dict = dict()
    mot_length = len(mask)
    logging.debug("AA motif table:\n%s", aa_df[(binom_prob <= args.p_value) & (fg_occ >= args.occurrences)].apply(lambda x: x.dropna().to_dict(), axis=1))
    for i in aa_df[(binom_prob <= args.p_value) & (fg_occ >= args.occurrences)].apply(lambda x: x.dropna().to_dict(), axis=1):
        if i:     
            for pos, aa in i.items():
                motif = mask[:]
                motifs_dict[pos] = aa
                index = pos + args.interval_length
                motif[index] = aa
                motifs.append([''.join(motif).strip('.'), binom_prob.at[aa, pos], fg_occ.at[aa, pos]])
    logging.debug("Motifs %s, mask %s", pd.DataFrame(motifs), mask)
    logging.debug("Motifs dict %s", motifs_dict)
    return motifs, motifs_dict


def get_df_mask(df, mask): 
    """
    Parameters
    ----------
    df : DataFrame
        DF with intervals
    mask : dict
        masks to apply

    Returns
    -------
    f_mask : Series
        combined mask

    """
    f_mask = True
    for i in [df[k] == v.upper() for k, v in mask.items()]:
        f_mask &= i
    return f_mask

 
def get_multimotifs(multi_list, single_dict):
    motifs_next = []
    for m_motif in multi_list:
        for s_motif, v in single_dict.items():
            if s_motif not in list(m_motif.keys()):
                motif = dict()
                motif.update({s_motif: v})
                motif.update(m_motif)
                if motif not in motifs_next:
                    motifs_next.append(motif)
    return motifs_next
   
def process_binom_stage(fg_intervals, bg_intervals, mask, results_saving_dir, args):
    motif_l = args.interval_length * 2 + 1
    logging.debug('FG %s mask is used\n%s', mask, fg_intervals[get_df_mask(fg_intervals, mask)] )
    fg_occ = utils.get_occurences(fg_intervals[get_df_mask(fg_intervals, mask)])
    if fg_occ[~fg_occ.isna()].sum().sum() > 0:
        logging.debug("Occurance fg matrix:\n%s", fg_occ)
        logging.debug('BG %s mask is used\n%s', mask, bg_intervals[get_df_mask(bg_intervals, mask)] )
        bg_occ = utils.get_occurences(bg_intervals[get_df_mask(bg_intervals, mask)])
        logging.debug("Occurance bg matrix:\n%s", bg_occ)
        binom_prob = get_binom_prob(fg_occ, bg_occ, args, results_saving_dir)
        motif = ['.'] * motif_l
        for pos, aa in mask.items():
            motif[pos + args.interval_length] = aa
        logging.debug('Prepared mask %s', motif)
        motifs, motif_dict = find_motifs(binom_prob, fg_occ, motif, args)
        return motifs, motif_dict
    else:
        return [], {}

def binomial_alg(fg_intervals, bg_intervals, args, results_saving_dir):
    final_motifs = []
    motifs, motifs_dict = process_binom_stage(fg_intervals, bg_intervals, {0: args.modification_site.lower()}, results_saving_dir, args)
    final_motifs.extend(motifs)
    multi_motifs = get_multimotifs([{k:v} for k, v in motifs_dict.items()], motifs_dict)
    logging.debug('Double motifs candidates %s', multi_motifs)
    while True:
        logging.debug('================================================NEW ITERATION===================================================================')
        tmp_mods_dict = dict()
        for motif in multi_motifs:
            new_motifs, new_motifs_dict = process_binom_stage(fg_intervals, bg_intervals, motif, results_saving_dir, args)
            final_motifs.extend(new_motifs)
            tmp_mods_dict.update(new_motifs_dict)
        if tmp_mods_dict:
            multi_motifs = get_multimotifs(tmp_mods_dict, motifs_dict.items())
        else:
            break
    print(pd.DataFrame(final_motifs))
    
