# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:36:07 2020

@author: Lab006
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


parser = argparse.ArgumentParser(description='PTM identification in modified peptides')
parser.add_argument('dataset', type=str, help ='Path to experimental dataset')
parser.add_argument('fasta',type=str, help='Path to FASTA')
parser.add_argument('--name_sample', type=str,default='sample_1', help='Name of examined sample',required=False)
parser.add_argument('--interval_length',type=int,default=6,help='Number of amino acids before & after modified amino acid (default=6)',required=False)
parser.add_argument('--modification',type=str,default='modification_1',help='Name of modification(ex.PHOSPHORYLATION)',required=False)
parser.add_argument('--modification_site',type=str,default='S',help='Modified amino acid (ex.S,T)')
parser.add_argument('--working_dir',type=str,default='.',help='Working dir for program (default=".")',required=False)
parser.add_argument('--algorithm',type=str,default='binom',help='Enter algorithm name: binom or chi2(default="binom")',required=False)
parser.add_argument('-v', '--verbosity', type=int, choices=range(3), default=1, help='Output verbosity',required=False)
args = parser.parse_args()


sample_saving_dir,results_saving_dir=utils.saving(args.working_dir,args.name_sample,args.interval_length,args.modification,args.modification_site,args.algorithm)    
#logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename=os.path.join(results_saving_dir,'mylog.log'))
levels = [logging.WARNING, logging.INFO, logging.DEBUG]
logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level=levels[args.verbosity])
logging.info(msg=u'Directories for result saving are created')


path_FASTA=args.fasta    


experimental_dir=args.dataset



#создаем таблицу пептидов на вход в программу
xl = pd.ExcelFile(experimental_dir)
df1 = xl.parse('Sheet1')
df1.columns=df1.values[0]
df=df1.drop(0)
Peptides=((df.loc[df['Amb?'] == 'UNIQUE']).reset_index()).loc[:,df.columns.intersection(['Peptide']) ]
indexes=[elem.replace('*','') for elem in Peptides['Peptide']]
Peptides['Mod_Peptide']=indexes  
Peptides.columns = ['Mod_Peptide','Peptide']
Peptides['Protein']=None
ind=Peptides['Peptide'].values
Peptides['index']=ind
Peptides=Peptides.set_index('index')

logging.debug(msg=u'Modified peptides dataframe is created')







a,b,c,d=utils.output(args.algorithm,args.working_dir,args.name_sample,Peptides,args.interval_length,args.modification_site,args.modification,path_FASTA)