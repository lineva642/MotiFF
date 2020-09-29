# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:35:59 2020

@author: Лиля
"""

import os
#from urllib.request import urlretrieve
#import pylab
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
#from pyteomics import fasta
#import gzip
#import seaborn as sns
#from scipy.stats import binom
import math
#import random
#import profile
from scipy.stats import chisquare

#import argparse
import logging
import utils

##АЛГОРИТМ_chi2

def  chi2_calc(occurrences,norm_expected,interval_length):
    p_value=pd.DataFrame(index=occurrences.index,columns=occurrences.columns)
    all_occ=(occurrences[interval_length]).sum()
    all_exp=(norm_expected[interval_length]).sum()
    for i in occurrences.columns:
        for j in occurrences.index:
            chisq,p_value[i][j]=chisquare([occurrences[i][j],all_occ-occurrences[i][j]],f_exp=[norm_expected[i][j],all_exp-norm_expected[i][j]])
    return p_value        
# считаем p-value
def p_value(occurrences,background_n,interval_length,results_saving_dir):
    
    #cчитаем количество интервалов в background
    all_back=(background_n[interval_length]).sum()
    #считаем количество интервалов в экспериментальном наборе
    all_exp=(occurrences[interval_length]).sum()
    #рассчитываем ожидаемую частоту встречаемости по распределению в background
    norm_expected=(background_n/all_back)*all_exp
    
    p_value=chi2_calc(occurrences,norm_expected,interval_length)
#    print(p_value)
                
    #результат нужно сохранить
    p_value.to_csv(os.path.join(results_saving_dir, 'p_value.csv'), sep='\t')

    logging.debug(msg=u'p-value matrix was created')
    return p_value

def primary_motifs(occurrences,background_n,p_value,acids,args,results_saving_dir):

    result=p_value[p_value<args.p_value/occurrences.size]*occurrences[occurrences>args.occurrences]
#    print(result)
    primary_motifs_number=[]
    primary_motifs_letter=[]
    for i in result.columns:
        for j in result.index:
            if (math.isnan(result[i][j])):
                continue
            else:

                if i<0:
                    n_motif=np.zeros(args.interval_length*2+1)
                    n_motif[i+args.interval_length]=acids.index(j)+1
                    motif=j+'.'*(-i-1)+args.modification_site.lower()
                    
                elif i>0:
                    n_motif=np.zeros(args.interval_length*2+1)
                    n_motif[i+args.interval_length]=acids.index(j)+1
                    motif=args.modification_site.lower()+'.'*(i-1)+j
                primary_motifs_letter.append(motif)
                primary_motifs_number.append(n_motif)
    vector=np.array(primary_motifs_number)    
    table=pd.DataFrame({'Number motif':primary_motifs_number,'Letter motif':primary_motifs_letter})      
    utils.saving_table(results_saving_dir,table,args.interval_length,'primary')
    logging.info(str(len(table['Number motif'].values))+u' primary (one acid length) motifs were identificated')
#    print(table)            
    return vector,table         



def counter(args, acid_location, acid_number, dataset_info, acids=utils.ACIDS_LIST):

    int_table, back_table , int_len, back_len = dataset_info 
    for i,ik in zip(acid_location, acid_number):
        int_table=int_table[(int_table[i-args.interval_length]==acids[ik-1])]
        back_table=back_table[(back_table[i-args.interval_length]==acids[ik-1])]
    observed=len(int_table)
    fasta = len(back_table)
    expected=(fasta/back_len)*int_len
    return observed,expected

    
def letter_motif(args,acid_location, acid_number,acids=utils.ACIDS_LIST):
    position=dict()
    for i,ik in zip(acid_location, acid_number):
        position[i]=acids[ik-1]
    position[args.interval_length]=(args.modification_site).lower()    

    keys=list(position.keys())
    keys.sort()
    motif=position[keys[0]]+'.'*(keys[1]-keys[0]-1)
    i=1
    while i<len(keys)-1:
        motif=''.join([motif,position[keys[i]],'.'*(keys[i+1]-keys[i]-1)])
        i+=1
    motif=''.join([motif,position[keys[len(keys)-1]]])
    return motif    
           
def chi2_motifs(args, acid_location, acid_number , dataset_info):
    

    observed,expected = counter(args, acid_location, acid_number, dataset_info, acids=utils.ACIDS_LIST)
    chisq,p_value=chisquare([observed,dataset_info[2]-observed],f_exp=[expected,dataset_info[2]-expected])

    return observed, expected, p_value

def acids_n_l(motif):
    
    acid_location = [elem for elem in np.nonzero(motif)[0]]
#                acid_number = [int(motif[np.nonzero(motif)[0][0]]),int(motif[np.nonzero(motif)[0][1]])]
    acid_number = [int(motif[elem]) for elem in np.nonzero(motif)[0]]
    
    return acid_location, acid_number

def multiplying_matrix(args, vector_1, vector_2):
    
    b=np.tensordot(vector_1,vector_2.T,axes=0)
    matrix=np.zeros((len(vector_1),len(vector_2),args.interval_length*2+1))

    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j]) 
    return matrix

def table_creator():
    
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])}) 
    return result

def motif_counter(args,motif,result,dataset_info):
    acid_location, acid_number = acids_n_l(motif)
    observed, expected, p_value=chi2_motifs(args, acid_location, acid_number , dataset_info)
    motif_l=letter_motif(args,acid_location, acid_number,acids=utils.ACIDS_LIST)
    result=result.append({'Letter motif':motif_l,'Number motif':motif,
                    'Observed':observed,'Expected':expected,
                                        'p-value':p_value},
                                        ignore_index=True)
    return result
  

def double_motifs(vector, dataset_info, results_saving_dir, args, acids=utils.ACIDS_LIST): 
    
    matrix = multiplying_matrix(args, vector, vector)
    print('double_motifs',vector)                   
    #создаем пустую табличку, в которую будем записывать результаты
    result = table_creator()
   
    for i in range(len(vector)):
        j=0
        while j<=i:
            elem=matrix[i,j]
            #нужны элементы матрицы с одними нулями
            if (elem.any())==False:
                motif=vector[i]+vector[j]
                result = motif_counter(args,motif,result,dataset_info)

            j+=1           
#    print(result)
    utils.saving_table(results_saving_dir,result,args.interval_length,'double')
    logging.info(str(len(result['Letter motif'].values))+u' double (two acid length) motifs were identificated')                                                                       
    return result

 
   
def n_vectors(args, previous_motifs_table):
    
    b=len(previous_motifs_table['Observed'].values)
    table=(previous_motifs_table.loc[previous_motifs_table['p-value']<args.p_value/b][previous_motifs_table['Observed']>=args.occurrences]).reset_index()
    del table['index']
    
    vector=np.zeros((len(table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        vector[i,:]=table['Number motif'].values[i]
        
    return vector

def n_motifs_result(args, n_vector, vector, matrix, result, dataset_info):
    for i in range(len(n_vector)):
        for j in range(len(vector)):
            elem=matrix[i,j]

            if (elem.any())==False:
                motif=vector[j]+n_vector[i]
                
              
                result = motif_counter(args,motif,result,dataset_info)
    return result            
    

def triple_motifs(vector,double_motifs,dataset_info,results_saving_dir, args, acids=utils.ACIDS_LIST):
       
    double_vector =  n_vectors(args, double_motifs) 
    matrix = multiplying_matrix(args, double_vector, vector)            
    result = table_creator()    
    result = n_motifs_result(args, double_vector, vector, matrix, result, dataset_info)
    utils.saving_table(results_saving_dir,result,args.interval_length,'triple')               
    logging.info(str(len(result['Letter motif'].values))+u' triple (three acid length) motifs were identificated')
        
    return result


#делаем прогу для 4хбуквенного мотива
def quadruple_motifs(vector, triple_motifs, dataset_info, results_saving_dir, args, acids=utils.ACIDS_LIST):
    
    triple_vector =  n_vectors(args, triple_motifs)
    matrix = multiplying_matrix(args, triple_vector, vector)    
    result = table_creator()                
    result = n_motifs_result(args, triple_vector, vector, matrix, result, dataset_info)
    utils.saving_table(results_saving_dir,result,args.interval_length,'quadruole')                    
    logging.info(str(len(result['Letter motif'].values))+u' quadruple (four acid length) motifs were identificated')
    
    return result


# In[117]: 


#улучшенная версия 2
    
def motifs(idPeptides, background, occurrences,background_n,p_value,args,results_saving_dir):
    vector,single=primary_motifs(occurrences,background_n,p_value,utils.ACIDS_LIST,args,results_saving_dir)
    intervals=(idPeptides['fasta_match']).sum()
    int_table = pd.DataFrame([list(i) for i in intervals], columns=range(-args.interval_length, args.interval_length + 1))
    back_table = pd.DataFrame([list(i) for i in background], columns=range(-args.interval_length, args.interval_length + 1))
    back_len=len(background)
    int_len=len(intervals)
    dataset_info = int_table, back_table , int_len, back_len    
    double=double_motifs(vector, dataset_info, results_saving_dir, args, acids=utils.ACIDS_LIST)
    if double is not None:
        triple=triple_motifs(vector,double,dataset_info,results_saving_dir,args, acids=utils.ACIDS_LIST)
        if triple is not None:
            quadruple=quadruple_motifs(vector, triple, dataset_info, results_saving_dir, args, acids=utils.ACIDS_LIST)
        else:
            quadruple=None
    else:
        triple=None   
    return single, double, triple, quadruple
