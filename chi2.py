# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:35:59 2020

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
    print(p_value)
                
    #результат нужно сохранить
    p_value.to_csv(os.path.join(results_saving_dir, 'p_value.csv'), sep='\t')

    logging.debug(msg=u'p-value matrix was created')
    return p_value

def primary_motifs(occurrences,background_n,p_value,acids,modification_site,interval_length,results_saving_dir):
#    p_value_selected=p_value[p_value<0.05/occurrences.size]
#    occurrences_selected=occurences[occurrences>10]
    result=p_value[p_value<0.05/occurrences.size]*occurrences[occurrences>10]
    print(result)
    primary_motifs_number=[]
    primary_motifs_letter=[]
    for i in result.columns:
        for j in result.index:
            if (math.isnan(result[i][j])):
                continue
            else:
#                print((i,j),result[i][j],type(result[i][j]))
                if i<0:
                    n_motif=np.zeros(interval_length*2+1)
                    n_motif[i+interval_length]=acids.index(j)+1
                    motif=j+'.'*(-i-1)+modification_site.lower()
                    
                elif i>0:
                    n_motif=np.zeros(interval_length*2+1)
                    n_motif[i+interval_length]=acids.index(j)+1
                    motif=modification_site.lower()+'.'*(i-1)+j
                primary_motifs_letter.append(motif)
                primary_motifs_number.append(n_motif)
    vector=np.array(primary_motifs_number)    
    table=pd.DataFrame({'Number motif':primary_motifs_number,'Letter motif':primary_motifs_letter})      
    utils.saving_table(results_saving_dir,table,interval_length,'primary')
    logging.info(str(len(table['Number motif'].values))+u' primary (one acid length) motifs were identificated')
    print(table)            
    return vector,table         


# In[114]:
def counter_doubles(dataset, i_1, i_2, k_1, k_2, acids=utils.ACIDS_LIST):
    count=0
    AA_1,AA_2=acids[i_1-1],acids[i_2-1]
    for elem in dataset:
        if (elem[k_1]==AA_1) and (elem[k_2]==AA_2):
            count+=1
    return count

def double_motifs(vector, occurrences, background_n, results_saving_dir, interval_length, modification_site, acids=utils.ACIDS_LIST):
    #выполняем тензорное умножение
    b = np.tensordot(vector, vector.T, axes=0)
    l=len(vector)    
    matrix=np.zeros((l, l, interval_length * 2 + 1))
    for i,ik in enumerate(b):    
        j= 0    
        while j<= i:        
            matrix[i,j,:] = np.diag(b[i,:,:,j])
            j+=1                
    
    #создаем пустую табличку, в которую будем записывать результаты
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})     
    back_len=len(background)
    int_len=len(intervals)
    for i in range(l):
        j=0
        while j<=i:
            elem=matrix[i,j]
            if (elem.any())==False:
                motif=vector[i]+vector[j]
                #положение аминокислот в интервале
                k_1,k_2=np.nonzero(motif)[0][0],np.nonzero(motif)[0][1]
                #номера аминокислот
                i_1,i_2=int(motif[k_1]),int(motif[k_2])
                observed=counter_doubles(intervals, i_1, i_2, k_1, k_2, acids=utils.ACIDS_LIST)
                fasta=counter_doubles(background, i_1, i_2, k_1, k_2, acids=utils.ACIDS_LIST)
#                observed=count(intervals,i_1,i_2,k_1,k_2,acids)
#                fasta=count(background,i_1,i_2,k_1,k_2,acids)
                expected=(fasta/back_len)*int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])
                position={k_1:acids[i_1-1],k_2:acids[i_2-1],interval_length:modification_site.lower()}
                keys=list(position.keys())
                keys.sort()
                motif_l=''.join([position[keys[0]],'.'*(keys[1]-keys[0]-1),position[keys[1]],'.'*(keys[2]-keys[1]-1),position[keys[2]]])
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True) 
            j+=1           
    logging.info(str(len(result['Letter motif'].values))+u' double (two acid length) motifs were identificated')                                                                       
    return result

def counter_triples(dataset, i_1, i_2, i_3, k_1, k_2, k_3, acids=utils.ACIDS_LIST):
    count=0
    AA_1,AA_2,AA_3=acids[i_1-1],acids[i_2-1],acids[i_3-1]
    for elem in dataset:
        if (elem[k_1]==AA_1) and (elem[k_2]==AA_2) and (elem[k_3]==AA_3):
            count+=1
    return count

def triple_motifs(primary_motifs,double_motifs,intervals,background,interval_length,results_saving_dir,modification_site, acids=utils.ACIDS_LIST):
    #составили всевозможные пары тройных мотивов
#    table=((double_motifs.where(double_motifs['Observed']>=10)).dropna()).reset_index()
    b=len(double_motifs['Observed'].values)
    table=(double_motifs.loc[double_motifs['p-value']<0.05/b][double_motifs['Observed']>=10]).reset_index()
    del table['index']

    double_vector=np.zeros((len(table['Number motif'].values),interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        double_vector[i,:]=table['Number motif'].values[i]
#    print(double_vector[0])

    vector=np.zeros((len(primary_motifs['Number motif'].values),interval_length*2+1)) 
    for i,ik in enumerate(primary_motifs['Number motif'].values):
        vector[i,:]=primary_motifs['Number motif'].values[i]
#    print(vector.shape)
    
    b=np.tensordot(double_vector,vector.T,axes=0)
    matrix=np.zeros((len(double_vector),len(vector),interval_length*2+1))
    back_len=len(background)
    int_len=len(intervals)
    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j])
            
#    print(matrix.shape)
#    print(matrix[0][0],vector[0],double_vector[0])
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})     
    for i in range(len(double_vector)):
        for j in range(len(vector)):
            elem=matrix[i,j]
#            print('elem',elem)
            if (elem.any())==False:
                motif=vector[j]+double_vector[i]
#                print('motif',motif)
#                print('vector',vector[j])
#                print('double_vector',double_vector[i])
                #положение аминокислот в интервале
                
                k_1,k_2,k_3=np.nonzero(motif)[0][0],np.nonzero(motif)[0][1],np.nonzero(motif)[0][2]
                #номера аминокислот
                i_1,i_2,i_3=int(motif[k_1]),int(motif[k_2]),int(motif[k_3])
#                print('!',k_1,k_2,k_3,i_1,i_2,i_3)
                observed=counter_triples(intervals, i_1, i_2, i_3, k_1, k_2, k_3, acids=utils.ACIDS_LIST)
                fasta=counter_triples(background, i_1, i_2, i_3, k_1, k_2, k_3, acids=utils.ACIDS_LIST)
                expected=(fasta / back_len) * int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])
                position={k_1:acids[i_1-1],k_2:acids[i_2-1],k_3:acids[i_3-1],interval_length:modification_site.lower()}
                keys=list(position.keys())
                keys.sort()
                motif_l=''.join([position[keys[0]],'.'*(keys[1]-keys[0]-1),position[keys[1]],'.'*(keys[2]-keys[1]-1),position[keys[2]],'.'*(keys[3]-keys[2]-1),position[keys[3]]])
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True)               
            
            
    logging.info(str(len(result['Letter motif'].values))+u' triple (three acid length) motifs were identificated')        
    return result

def counter_quadruples(dataset, i_1, i_2, i_3, i_4, k_1, k_2, k_3, k_4, acids=utils.ACIDS_LIST):
    count=0
    AA_1,AA_2,AA_3,AA_4=acids[i_1-1],acids[i_2-1],acids[i_3-1],acids[i_4-1]
    for elem in dataset:
        if (elem[k_1]==AA_1) and (elem[k_2]==AA_2) and (elem[k_3]==AA_3) and (elem[k_4]==AA_4):
            count+=1
    return count
#делаем прогу для 4хбуквенного мотива
def quadruple_motifs(primary_motif, triple_motifs, intervals, background, interval_length, results_saving_dir, modification_site, acids=utils.ACIDS_LIST):
    
    b=len(triple_motifs['Observed'].values)
    table=(triple_motifs.loc[triple_motifs['p-value']<0.05/b][triple_motifs['Observed']>=10]).reset_index()
    del table['index']
#    print(table)
    
    triple_vector=np.zeros((len(table['Number motif'].values),interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        triple_vector[i,:]=table['Number motif'].values[i]
#    print(triple_vector[0])

    vector=np.zeros((len(primary_motif['Number motif'].values),interval_length*2+1)) 
    for i,ik in enumerate(primary_motif['Number motif'].values):
        vector[i,:]=primary_motif['Number motif'].values[i]
#    print(vector.shape)
    
    b=np.tensordot(triple_vector,vector.T,axes=0)
    matrix=np.zeros((len(triple_vector),len(vector),interval_length*2+1))
    back_len=len(background)
    int_len=len(intervals)
    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j])
            
#    print(matrix.shape)
#    print(matrix[0][0],vector[0],double_vector[0])
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})
    for i in range(len(triple_vector)):
        for j in range(len(vector)):
            elem=matrix[i,j]
#            print('elem',elem)
            if (elem.any())==False:
                motif=vector[j]+triple_vector[i]
#                print('motif',motif)
#                print('vector',vector[j])
#                print('double_vector',double_vector[i])
                #положение аминокислот в интервале
                
                k_1,k_2,k_3,k_4=np.nonzero(motif)[0][0],np.nonzero(motif)[0][1],np.nonzero(motif)[0][2],np.nonzero(motif)[0][3]
                #номера аминокислот
                i_1,i_2,i_3,i_4=int(motif[k_1]),int(motif[k_2]),int(motif[k_3]),int(motif[k_4])
#                print('!',k_1,k_2,k_3,k_4,i_1,i_2,i_3,i_4)
                observed=counter_quadruples(intervals, i_1, i_2, i_3, i_4, k_1, k_2, k_3, k_4, acids=utils.ACIDS_LIST)
                fasta=counter_quadruples(background, i_1, i_2, i_3, i_4, k_1, k_2, k_3, k_4, acids=utils.ACIDS_LIST)
                expected=(fasta/back_len)*int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])
                position={k_1:acids[i_1-1],k_2:acids[i_2-1],k_3:acids[i_3-1],k_4:acids[i_4-1],interval_length:modification_site.lower()}
                keys=list(position.keys())
                keys.sort()
                motif_l=''.join([position[keys[0]],'.'*(keys[1]-keys[0]-1),position[keys[1]],'.'*(keys[2]-keys[1]-1),position[keys[2]],'.'*(keys[3]-keys[2]-1),position[keys[3]],
                                 '.'*(keys[4]-keys[3]-1),position[keys[4]]])
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True)                    
    logging.info(str(len(result['Letter motif'].values))+u' quadruple (four acid length) motifs were identificated')
    return result


# In[117]: 


#улучшенная версия 2
#def motifs(chi2_selection,interval_length,modification_site,background,intervals,results_saving_dir, acids=utils.ACIDS_LIST):
def motifs(occurrences,background_n,p_value,args,results_saving_dir):
#    vector,single=primary_motifs(interval_length, chi2_selection, results_saving_dir, modification_site, acids=utils.ACIDS_LIST)
#    p_value_selected=primary_motifs(occurrences,background_n,p_value)
    vector,table=primary_motifs(occurrences,background_n,p_value,utils.ACIDS_LIST,args.modification_site,args.interval_length,results_saving_dir)
#    double=double_motifs(vector, intervals, background, results_saving_dir, interval_length, modification_site, acids=utils.ACIDS_LIST)
#    if double is not None:
#        triple=triple_motifs(single, double, intervals, background, interval_length, results_saving_dir, modification_site, acids=utils.ACIDS_LIST)
#        if triple is not None:
#            quadruple=quadruple_motifs(single, triple, intervals, background, interval_length, results_saving_dir, modification_site, acids=utils.ACIDS_LIST)
#        else:
#            quadruple=None
#    else:
#        triple=None   
#    return single, double, triple, quadruple
    return vector,table
#    return single,double