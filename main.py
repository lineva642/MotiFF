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
parser.add_argument('fasta', type=str, help='Path to FASTA')
parser.add_argument('--name_sample', type=str,default='sample_1', help='Name of examined sample', required=False)
parser.add_argument('--interval_length', type=int, default=4, 
                    help='Number of amino acids before & after modified amino acid (default=6)', required=False)
parser.add_argument('--modification', type=str, default='modification_1', 
                    help='Name of modification(ex.PHOSPHORYLATION)', required=False)
parser.add_argument('--modification_site', type=str, default='S', help='Modified amino acid (ex.S,T)')
parser.add_argument('--working_dir', type=str, default='.', help='Working dir for program (default=".")', required=False)
parser.add_argument('--algorithm', type=str, default='chi2',help='Enter algorithm name: binom or chi2(default="chi2")', required=False)
parser.add_argument('--p_value', type=float, default=0.000005, help='Enter p_value(default=0.05)', required=False)
parser.add_argument('--occurrences', type=int, default=20,
                    help='Enter number of motif occurrences in experimental dataset (default=10)', required=False)
parser.add_argument('-v', '--verbosity', type=int, choices=range(3), default=1, help='Output verbosity', required=False)
args = parser.parse_args()


levels = [logging.WARNING, logging.INFO, logging.DEBUG]
logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level=levels[args.verbosity])
logging.info(msg=u'Directories for result saving are created')

experimental_dir = args.dataset

logging.debug(msg=u'Modified peptides dataframe is created')


a,b,c,d=utils.output(args)