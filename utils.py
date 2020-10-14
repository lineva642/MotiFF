# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:37:30 2020

@author: Lab006

"""

ACIDS_LIST = ['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D','-']
import os
from urllib.request import urlretrieve
import pylab
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyteomics import fasta
import gzip
#import seaborn as sns
# from scipy.stats import binom
import math
import random
import profile
from scipy.stats import chisquare

import argparse
import logging
import chi2
import binomial
from collections import defaultdict, Counter
import re


def saving(args):

    sample_saving_dir = os.path.join(args.working_dir, args.name_sample)
    results_saving_dir = os.path.join(sample_saving_dir, args.modification + '_' +
                                      args.modification_site + str(args.interval_length) + '_' + args.algorithm)

    if not (os.path.exists(sample_saving_dir)):
        os.mkdir(sample_saving_dir)
        if not (os.path.exists(results_saving_dir)):
            os.mkdir(results_saving_dir)
    else:
        if not (os.path.exists(results_saving_dir)):
            os.mkdir(results_saving_dir)
#    logging.basicConfig(level = logging.DEBUG,filename=os.path.join(results_saving_dir,'mylog.log'))        
#    logging.info(u'Directories for result saving are created')       
    return  sample_saving_dir, results_saving_dir


def fasta_match(row, bg_fasta, interval_length, modification_site):
    intervals = []
    # k = 0
    # print(row['Peptide'])
    for name, seq in bg_fasta.items():
        i = 0
        start = seq[i:].find(row['Peptide'].replace('*', ''))
        i = start
        while start >= 0:
            for asterisks, modif in enumerate(re.finditer('\*', row['Peptide']), 1):
                interval_start = i + modif.span()[0] - interval_length - asterisks
                interval_end = interval_start + 2 * interval_length + 1
                interval = seq[interval_start: interval_end]
                if interval[interval_length] == modification_site:
                    intervals.append(interval)
                else:
                    logging.info('Peptide has another modification site %s', row['Peptide'])
            start = seq[i + 1:].find(row['Peptide'].replace('*', ''))
            i += start + 1
    # k += 1        
    return intervals


def background_maker(args):
#    print('Making background DB')
    #хотим сделать background из идентифицированных белков
    bg_fasta = dict()
    # bg = defaultdict()
    background = set()
    with fasta.read(args.fasta) as f:
        for name, sequence in f:
            name_id = name.split('|')[1]
            extended_seq = ''.join(['-' * args.interval_length, sequence, '-' * args.interval_length])
            bg_fasta[name_id] = extended_seq
            mod_aa_indexes = re.finditer(args.modification_site, extended_seq)
            bg_intervals = [extended_seq[i.span()[0] - args.interval_length: i.span()[0] + args.interval_length + 1] for i in mod_aa_indexes]
            # bg[name_id] = bg_intervals
            background.update(bg_intervals)

    logging.info(u'Set of %s background intervals is created', len(background))
    logging.debug(u'Background DB is ready')    
    return background, bg_fasta   


#функции для валидации

def aa_counter(col):
    return pd.Series(Counter(col))


def get_occurences(intervals_list, interval_length, saving_file, acids=ACIDS_LIST):
    logging.debug('Intervals list length :\n%s', len(intervals_list))
    df = pd.DataFrame([list(i) for i in intervals_list], columns=range(-interval_length, interval_length + 1))
    occ = pd.DataFrame(index=acids)
    occ[df.columns] = df.apply(aa_counter, axis=0)
    logging.debug("Occurance matrix (%s):\n%s", saving_file, occ)
    occ.to_csv(saving_file, sep='\t')
    return occ


def saving_table(results_saving_dir, result, interval_length, name):
    path=os.path.join(results_saving_dir, 'table' + str(interval_length) + '_' + name + '.csv')
    result.to_csv(path)   
 
    
def peptides_table(args, sample_saving_dir, bg_fasta):
    # if os.path.exists(os.path.join(sample_saving_dir, 'peptide_identification.csv')):
        
    #     idPeptides = pd.read_csv(os.path.join(sample_saving_dir, 'peptide_identification.csv'), sep=',', usecols=[1,2])
    #     for i in range(len(idPeptides['fasta_match'].index)):
    #         idPeptides['fasta_match'][i]=(((idPeptides['fasta_match'][i].replace('[','')).replace(']','')).replace("'",'')).split(',')
    # else:    
        # peptides = []
        # with open(args.dataset) as f:
        #     for line in f:
        #         peptides.append(line[:-1])    
        # Peptides = pd.DataFrame({'Peptide':peptides})
    Peptides = pd.read_table(args.dataset, names=['Peptide'])
    logging.info('Peptide table contains %d peptides', Peptides.shape[0])
    logging.debug('Initial Peptides df:\n%s', Peptides.head())
    Peptides['fasta_match'] = Peptides.apply(fasta_match, args=[bg_fasta, args.interval_length, args.modification_site], axis=1)
    # print(background)np
    Peptides['unique'] = Peptides.apply(lambda x: True if len(x['fasta_match']) == 1 else False, axis=1)
    # print(Peptides)
    idPeptides = Peptides[Peptides['unique'] == True]
    logging.info('Found %d unique peptides motifs', len(idPeptides))
    idPeptides.to_csv(os.path.join(sample_saving_dir, 'peptide_identification.csv'), mode='w')
    logging.debug('Prepared Peptides df:\n%s', idPeptides) 
    return idPeptides    


def output(args):
    
    sample_saving_dir, results_saving_dir = saving(args)
    background, bg_fasta = background_maker(args)
    idPeptides = peptides_table(args, sample_saving_dir, bg_fasta)
    
    occurrences = get_occurences((idPeptides['fasta_match']).sum(), args.interval_length,
                             os.path.join(results_saving_dir, 'occurences.csv'), 
                             acids=ACIDS_LIST)
    background_n = get_occurences(background, args.interval_length, 'background.csv')
    logging.debug("Final idPeptides table:\n%s", idPeptides.head())

    if args.algorithm == "binom":
        logging.debug('Binomial algorithm is used.')        
        P_binomial = binomial.P_counter_bi(occurrences, background_n, args, results_saving_dir)
        single, double, triple, quadruple = binomial.motifs_bi(P_binomial, occurrences, idPeptides, background, args, results_saving_dir)

        logging.info(msg='Program was finished successfully') 
        return single,double,triple,quadruple
    else:
        
#        occurrences = get_occurences( (idPeptides['fasta_match']).sum(), args.interval_length,
#                                     os.path.join(results_saving_dir, 'occurences.csv'), 
#                                     acids=ACIDS_LIST)
#        background_n = get_occurences(background, args.interval_length, 'background.csv')

        p_value=chi2.p_value(occurrences,background_n,args.interval_length,results_saving_dir)

        single, double, triple, quadruple=chi2.motifs(idPeptides, background, occurrences,background_n,p_value,args,results_saving_dir)
        logging.info(msg='Program was finished successfully') 
#        return chi2_results,chi2_selection,intervals,background
        return single, double, triple, quadruple

