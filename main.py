# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:36:07 2020

@author: Lab006
"""

import pandas as pd
import argparse
import logging
import utils, binomial, chi2


parser = argparse.ArgumentParser(description='PTM identification in modified peptides')
parser.add_argument('dataset', type=str, help ='Path to experimental dataset')
parser.add_argument('fasta', type=str, help='Path to FASTA')
parser.add_argument('--name_sample', type=str,default='sample_1', help='Name of examined sample', required=False)
parser.add_argument('--interval_length', type=int, default=6, 
                    help='Number of amino acids before & after modified amino acid (default=6)', required=False)
parser.add_argument('--modification', type=str, default='modification_1', 
                    help='Name of modification(ex.PHOSPHORYLATION)', required=False)
parser.add_argument('--modification_site', type=str, default='S', help='Modified amino acid (ex.S,T)')
parser.add_argument('--working_dir', type=str, default='.', help='Working dir for program (default=".")', required=False)
parser.add_argument('--algorithm', type=str, default='chi2',help='Enter algorithm name: binom or chi2(default="chi2")', required=False)
parser.add_argument('--p_value', type=float, default=0.000005, help='Enter p_value(default=0.000005)', required=False)
parser.add_argument('--occurrences', type=int, default=20,
                    help='Enter number of motif occurrences in experimental dataset (default=20)', required=False)
parser.add_argument('-v', '--verbosity', type=int, choices=range(3), default=1, help='Output verbosity', required=False)
args = parser.parse_args()


levels = [logging.WARNING, logging.INFO, logging.DEBUG]
logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level=levels[args.verbosity])
logging.info(msg=u'Directories for result saving are created')

experimental_dir = args.dataset

logging.debug(msg=u'Modified peptides dataframe is created')

def output(args):
    
    sample_saving_dir, results_saving_dir = utils.saving(args)
    bg_intervals, bg_fasta = utils.background_maker(args)
    idPeptides = utils.peptides_table(args, sample_saving_dir, bg_fasta)
    # print("idPeptides['fasta_match'].sum()",idPeptides['fasta_match'],idPeptides['fasta_match'].sum())
    fg_intervals = pd.DataFrame([list(i) for i in idPeptides['fasta_match'].sum()], 
                                columns=range(-args.interval_length, args.interval_length + 1))
 
    logging.debug("Final idPeptides table:\n%s", idPeptides.head())

    if args.algorithm == "binom":
        logging.debug('Binomial algorithm is used.')  
        result = binomial.binomial_alg(fg_intervals, bg_intervals, args)
        result.to_csv('Binom_motifs.csv', sep='\t')
        
    else:
        dataset_info = fg_intervals, bg_intervals
        fg_occ = utils.get_occurences(fg_intervals, acids=utils.ACIDS_LIST)
        bg_occ = utils.get_occurences(bg_intervals, acids=utils.ACIDS_LIST)
        chi2.chi2_alg(fg_occ, bg_occ, dataset_info, args, results_saving_dir)

output(args)        
logging.info(msg='Program was finished successfully') 