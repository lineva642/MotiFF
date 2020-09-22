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


# In[114]:
#def counter(dataset, indexes, acids=utils.ACIDS_LIST):
#    count=0
#    for elem in dataset:
#        for index in indexes:
#            if elem[index]==acids[indexes[index]-1]:
#                x=True
#            else:
#                x=False
#                break    
#        if x==True:
#            count+=1
##    print(count)    
#    return count
    
#def counter(args, dataset, indexes, acids=utils.ACIDS_LIST):
#    table = pd.DataFrame([list(i) for i in dataset], columns=range(-args.interval_length, args.interval_length + 1))
#    for index in indexes:
#        table=table[table[index-args.interval_length]==acids[indexes[index]-1]]
#    count=len(table)
##    print('!',count)
#    return count

def counter_double(args, dataset, i_1, i_2, k_1, k_2, acids=utils.ACIDS_LIST):
    table = pd.DataFrame([list(i) for i in dataset], columns=range(-args.interval_length, args.interval_length + 1))
    table=table[(table[k_1-args.interval_length]==acids[i_1-1]) & (table[k_2-args.interval_length]==acids[i_2-1])]
#    for index in indexes:
#        table=table[table[index-args.interval_length]==acids[indexes[index]-1]]
    count=len(table)
#    print('!',count)
    return count
    
def letter_motif(args,indexes,acids=utils.ACIDS_LIST):
    position=dict()
    for elem in indexes:
        position[elem]=acids[indexes[elem]-1]
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
           
                

def double_motifs(vector, idPeptides, background, results_saving_dir, args, acids=utils.ACIDS_LIST):
    #выполняем тензорное умножение
    b = np.tensordot(vector, vector.T, axes=0)
    l=len(vector)    
    matrix=np.zeros((l, l, args.interval_length * 2 + 1))
    for i,ik in enumerate(b):    
        j= 0    
        while j<= i:        
            matrix[i,j,:] = np.diag(b[i,:,:,j])
            j+=1 

                   
    #создаем пустую табличку, в которую будем записывать результаты
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})
    intervals=(idPeptides['fasta_match']).sum()
    
    int_len= len(intervals)
    back_len= len(background)
    
    for i in range(l):
        j=0
        while j<=i:
            elem=matrix[i,j]
            #нужны элементы матрицы с одними нулями
            if (elem.any())==False:
                motif=vector[i]+vector[j]
                 #положение аминокислот в интервале
                k_1,k_2=np.nonzero(motif)[0][0],np.nonzero(motif)[0][1]
                 #номера аминокислот
                i_1,i_2=int(motif[k_1]),int(motif[k_2])
                
                indexes=dict(zip([np.nonzero(motif)[0][0],np.nonzero(motif)[0][1]],[int(motif[np.nonzero(motif)[0][0]]),int(motif[np.nonzero(motif)[0][1]])]))

                
#                observed=counter(args, intervals, indexes, acids=utils.ACIDS_LIST)
#                fasta=counter(args, background, indexes, acids=utils.ACIDS_LIST)
                observed=counter_double(args, intervals, i_1, i_2, k_1, k_2, acids=utils.ACIDS_LIST)
                fasta=counter_double(args, background, i_1, i_2, k_1, k_2, acids=utils.ACIDS_LIST)

                expected=(fasta/back_len)*int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])

                motif_l=letter_motif(args,indexes,acids=utils.ACIDS_LIST)
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True)
            j+=1           
    print(result)
    utils.saving_table(results_saving_dir,result,args.interval_length,'double')
    logging.info(str(len(result['Letter motif'].values))+u' double (two acid length) motifs were identificated')                                                                       
    return result

def triple_motifs(primary_motifs,double_motifs,idPeptides,background,results_saving_dir, args, acids=utils.ACIDS_LIST):
    #составили всевозможные пары тройных мотивов

    b=len(double_motifs['Observed'].values)
    table=(double_motifs.loc[double_motifs['p-value']<args.p_value/b][double_motifs['Observed']>=args.occurrences]).reset_index()
    del table['index']


    double_vector=np.zeros((len(table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        double_vector[i,:]=table['Number motif'].values[i]


    vector=np.zeros((len(primary_motifs['Number motif'].values),args.interval_length*2+1)) 
    for i,ik in enumerate(primary_motifs['Number motif'].values):
        vector[i,:]=primary_motifs['Number motif'].values[i]

    
    b=np.tensordot(double_vector,vector.T,axes=0)
    matrix=np.zeros((len(double_vector),len(vector),args.interval_length*2+1))
    intervals=(idPeptides['fasta_match']).sum()
    back_len=len(background)
    int_len=len(intervals)
    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j])
            
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})     
    for i in range(len(double_vector)):
        for j in range(len(vector)):
            elem=matrix[i,j]

            if (elem.any())==False:
                motif=vector[j]+double_vector[i]

                indexes=dict(zip([np.nonzero(motif)[0][0],np.nonzero(motif)[0][1],np.nonzero(motif)[0][2]],[int(motif[np.nonzero(motif)[0][0]]),int(motif[np.nonzero(motif)[0][1]]),int(motif[np.nonzero(motif)[0][2]])]))
                observed=counter(args, intervals, indexes, acids=utils.ACIDS_LIST)
                fasta=counter(args, background, indexes, acids=utils.ACIDS_LIST)
                expected=(fasta / back_len) * int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])

                motif_l=letter_motif(args,indexes,acids=utils.ACIDS_LIST)
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True)               
                
    utils.saving_table(results_saving_dir,result,args.interval_length,'triple')               
    logging.info(str(len(result['Letter motif'].values))+u' triple (three acid length) motifs were identificated')        
    return result

#делаем прогу для 4хбуквенного мотива
def quadruple_motifs(primary_motif, triple_motifs, idPeptides, background, results_saving_dir, args, acids=utils.ACIDS_LIST):
    
    b=len(triple_motifs['Observed'].values)
    table=(triple_motifs.loc[triple_motifs['p-value']<args.p_value/b][triple_motifs['Observed']>=args.occurrences]).reset_index()
    del table['index']

    
    triple_vector=np.zeros((len(table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        triple_vector[i,:]=table['Number motif'].values[i]


    vector=np.zeros((len(primary_motif['Number motif'].values),args.interval_length*2+1)) 
    for i,ik in enumerate(primary_motif['Number motif'].values):
        vector[i,:]=primary_motif['Number motif'].values[i]

    
    b=np.tensordot(triple_vector,vector.T,axes=0)
    matrix=np.zeros((len(triple_vector),len(vector),args.interval_length*2+1))
    
    intervals=(idPeptides['fasta_match']).sum()
    back_len=len(background)
    int_len=len(intervals)
    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j])
            

    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})
    for i in range(len(triple_vector)):
        for j in range(len(vector)):
            elem=matrix[i,j]

            if (elem.any())==False:
                motif=vector[j]+triple_vector[i]
                indexes=dict(zip([np.nonzero(motif)[0][0],np.nonzero(motif)[0][1],np.nonzero(motif)[0][2],np.nonzero(motif)[0][3]],[int(motif[np.nonzero(motif)[0][0]]),int(motif[np.nonzero(motif)[0][1]]),int(motif[np.nonzero(motif)[0][2]]),int(motif[np.nonzero(motif)[0][3]])]))

                observed=counter(args, intervals, indexes, acids=utils.ACIDS_LIST)
                fasta=counter(args, background, indexes, acids=utils.ACIDS_LIST)
                expected=(fasta/back_len)*int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])

                motif_l=letter_motif(args,indexes,acids=utils.ACIDS_LIST)
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True)
    print(result)
    utils.saving_table(results_saving_dir,result,args.interval_length,'quadruole')                    
    logging.info(str(len(result['Letter motif'].values))+u' quadruple (four acid length) motifs were identificated')
    return result


# In[117]: 


#улучшенная версия 2
    
def motifs(idPeptides, background, occurrences,background_n,p_value,args,results_saving_dir):
    vector,single=primary_motifs(occurrences,background_n,p_value,utils.ACIDS_LIST,args,results_saving_dir)
    double=double_motifs(vector, idPeptides, background, results_saving_dir, args, acids=utils.ACIDS_LIST)
    if double is not None:
        triple=triple_motifs(single,double,idPeptides,background,results_saving_dir,args, acids=utils.ACIDS_LIST)
        if triple is not None:
            quadruple=quadruple_motifs(single, triple, idPeptides, background, results_saving_dir, args, acids=utils.ACIDS_LIST)
        else:
            quadruple=None
    else:
        triple=None   
    return single, double, triple, quadruple
