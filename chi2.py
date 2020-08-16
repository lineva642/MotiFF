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
#считаем ожидаемое распределение аминокислот
def expected_destribution(background_n,occurrences,interval_length):

    all_back=0
    for k in range(0,21):
        all_back+=background_n[k][0]

    all_exp=0
    for k in range(0,21):
        all_exp+=occurrences[k][0]

    norm_FASTA=[]
    for i in range(0,21):
        a=[0]*(interval_length*2+1)
        norm_FASTA.append(a)
    for i in range(0,21):
        for k in range(0,interval_length*2+1):
            norm_FASTA[i][k]=(background_n[i][k])/all_back

    expected_FASTA=[]
    for i in range(0,21):
        a=[0]*(interval_length*2+1)
        expected_FASTA.append(a)
    for i in range(0,21):
        for k in range(0,interval_length*2+1):
            expected_FASTA[i][k]=(norm_FASTA[i][k])*all_exp
    logging.debug(u'Expected destrbution was counted')        
    return expected_FASTA,all_exp


# считаем p-value
def chi2_result(occurrences_FASTA,expected_FASTA,all_exp,interval_length,results_saving_dir):
    chi2_results=[]
    for i in range(0,21):
        a=[0]*(interval_length*2+1)
        chi2_results.append(a)
    for i in range(0,21):
        for k in range(0,interval_length*2+1):
            if k!=interval_length:
                chisq,p_value=chisquare([occurrences_FASTA[i][k],all_exp-occurrences_FASTA[i][k]],f_exp=[expected_FASTA[i][k],all_exp-expected_FASTA[i][k]])
                chi2_results[i][k]=p_value
    #результат нужно сохранить
    path=os.path.join(results_saving_dir,'p_value.txt')
    saving=open(path,'w')
    for i in range(0,21):
        for k in range(interval_length*2+1):
            saving.write(str(chi2_results[i][k])+' ')
        saving.write('\n')
    saving.close()
    logging.debug(msg=u'Chi2 values matrix was created')
    return chi2_results            


# In[107]:


def p_value_selection(chi2_results,occurrences_FASTA,interval_length):
    chi2_selection=[]
    for i in range(0,21):
        a=[0]*(2*interval_length+1)
        chi2_selection.append(a)
    for i in range(0,21):
        for k in range(0,2*interval_length+1):
            if (chi2_results[i][k]<0.05/(21*(2*interval_length+1))) and (k!=interval_length) and (occurrences_FASTA[i][k]>10):
                chi2_selection[i][k]=1
    logging.debug(u'Selection of matrix elements was successful')            
    return chi2_selection            


# In[108]:


#напишем одну функцию для подсчета p-value


# In[109]:


def p_value(background_n,occurrences,interval_length,modification_site,acids,results_saving_dir):
    expected_FASTA,all_exp=expected_destribution(background_n,occurrences,interval_length)
    chi2_results=chi2_result(occurrences,expected_FASTA,all_exp,interval_length,results_saving_dir)
    chi2_selection=p_value_selection(chi2_results,occurrences,interval_length)
#    heatmap_visualization(chi2_selection,acids,interval_length,modification_site,'Отбор по p-value',results_saving_dir,name='p_value_selection')
    logging.info(u'Chi2 values were counted and selection with p=0.05/Bonferroni and occurrences>10 correction was performed')
    return expected_FASTA,chi2_results,chi2_selection

#функция для подсчета комлексных мотивов




def primary_motifs(acids,interval_length,chi2_selection,results_saving_dir,modification_site):
    primary_motifs=[]
    primary_location=[]
    for i in range(0,21):
        for k in range(0,2*interval_length+1):
            if chi2_selection[i][k]==1:
                primary_motifs.append(acids[i])
                primary_location.append((i,k))            
    primary_motif=pd.DataFrame({'Acid':primary_motifs,'Location':primary_location})
    
    motifs=[]
    for i in range(len(primary_motif['Location'].values)):
#         print(i)
        motif=dict()
        pair_1_x,pair_1_y=primary_motif['Location'].values[i]
        motif[pair_1_y]=acids[pair_1_x]
        motif[interval_length]=modification_site.lower()
        list_keys = list(motif.keys())
        list_keys.sort()
        word_motif=motif[list_keys[0]]+'.'*(list_keys[1]-list_keys[0]-1)+motif[list_keys[1]]
        motifs.append(word_motif)
        
    primary_motif['Motifs']=motifs
    utils.saving_table(results_saving_dir,primary_motif,interval_length,'primary')
#    print('Primary motifs are ready!')
    logging.info(str(len(primary_motif['Motifs'].values))+u' primary (one acid length) motifs were identificated')
    return primary_motif


# In[114]:


def double_motifs(primary_motif,acids,intervals,background,results_saving_dir,interval_length,modification_site):
    primary_motif_copy=primary_motif.copy()
    second_motifs_location=[]
    for i in range(len(primary_motif_copy)):
        for k in range(i+1,len(primary_motif_copy)):
            second_motifs_location.append((primary_motif_copy['Location'].values[i],primary_motif_copy['Location'].values[k]))
    
    second_motifs=pd.DataFrame({'Location':second_motifs_location})
    
    #для каждой пары аминокислот должны рассчитать observed и expected
    occurrences=[]
    expactations=[]
    p_values=[]
    for i in range(len(second_motifs['Location'].values)):
#         print(i)
        first_AA,second_AA=second_motifs['Location'].values[i]
        first_AA_x,first_AA_y=first_AA
        second_AA_x,second_AA_y=second_AA
        observed=0
        for elem in intervals:
            if (elem[first_AA_y]==acids[first_AA_x]) and (elem[second_AA_y]==acids[second_AA_x]):
                observed+=1
        occurrences.append(observed)  


        fasta=0        
        for elem in background:
            if (elem[first_AA_y]==acids[first_AA_x]) and (elem[second_AA_y]==acids[second_AA_x]):
                fasta+=1

        expected=(fasta/len(background))*len(intervals)
        expactations.append(expected)

        if observed>10:
            chisq,p_value=chisquare([observed,len(intervals)-observed],f_exp=[expected,len(intervals)-expected])
            p_values.append(p_value)
        else:
            p_values.append('-')
            
    second_motifs['Observed']=occurrences
    second_motifs['Expected']=expactations
    second_motifs['p-value']=p_values
    
    second_motifs_selection=second_motifs.where(second_motifs['p-value']!='-').dropna()
    
    count=[]
    for i in range(len(second_motifs_selection['p-value'].values)):
        if (second_motifs_selection['p-value'].values[i])<0.05/len(second_motifs['Location'].values):
            count.append(1)
        else:
            count.append(0)
    second_motifs_selection['Count']=count
    
    second_motifs_selection_p=second_motifs_selection.where(second_motifs_selection['Count']==1).dropna()
    second_motifs_selection_p_copy=second_motifs_selection_p.copy()
    del second_motifs_selection_p_copy['Count']
    
    motifs=[]
    for i in range(len(second_motifs_selection_p_copy['Location'].values)):
#         print(i)
        motif=dict()
        first_AA,second_AA=second_motifs_selection_p_copy['Location'].values[i]
        first_AA_x,first_AA_y=first_AA
        second_AA_x,second_AA_y=second_AA
        motif[first_AA_y]=acids[first_AA_x]
        motif[second_AA_y]=acids[second_AA_x]
        motif[interval_length]=modification_site.lower()
        list_keys = list(motif.keys())
        list_keys.sort()
        word_motif=motif[list_keys[0]]+'.'*(list_keys[1]-list_keys[0]-1)+motif[list_keys[1]]+'.'*(list_keys[2]-list_keys[1]-1)+motif[list_keys[2]]
        motifs.append(word_motif)
        
    second_motifs_selection_p_copy['Motifs']=motifs
    second_motifs_selection_p_copy=second_motifs_selection_p_copy.reset_index()    
    del second_motifs_selection_p_copy['index']
    
    utils.saving_table(results_saving_dir,second_motifs_selection_p_copy,interval_length,'double')
#    print('Double motifs are ready!')
    logging.info(str(len(second_motifs_selection_p_copy['Motifs'].values))+u' double (two acid length) motifs were identificated')
    
    return second_motifs_selection_p_copy


# In[115]:


def triple_motifs(primary_motif,second_motifs,acids,intervals,background,interval_length,results_saving_dir,modification_site):
    location=[]
    for i in range(len(second_motifs['Location'].values)):
        pair_1,pair_2=second_motifs['Location'].values[i]
        for k in range(len(primary_motif['Location'].values)):
            pair_3=primary_motif['Location'].values[k]
            location.append((pair_1,pair_2,pair_3))
            
    triple_motifs=pd.DataFrame({'Location':location})
    
    count=[]
    for i in range(len(triple_motifs['Location'].values)):
        pair_1,pair_2,pair_3=triple_motifs['Location'].values[i]
        pair_1_x,pair_1_y=pair_1
        pair_2_x,pair_2_y=pair_2
        pair_3_x,pair_3_y=pair_3
        if (pair_1_y!=pair_2_y) and (pair_1_y!=pair_3_y) and (pair_3_y!=pair_2_y):
            count.append(1)
        else:
            count.append(0)
    triple_motifs['Count']=count
    triple_motifs=(triple_motifs.where(triple_motifs['Count']==1).dropna()).reset_index()
    del triple_motifs['Count']
    del triple_motifs['index']
    
    #для каждой тройки аминокислот должны рассчитать observed и expected
    occurrences=[]
    expactations=[]
    p_values=[]
    for i in range(len(triple_motifs['Location'].values)):
#         print(i)
        pair_1,pair_2,pair_3=triple_motifs['Location'].values[i]
        pair_1_x,pair_1_y=pair_1
        pair_2_x,pair_2_y=pair_2
        pair_3_x,pair_3_y=pair_3
        observed=0
        for elem in intervals:
            if (elem[pair_1_y]==acids[pair_1_x]) and (elem[pair_2_y]==acids[pair_2_x]) and (elem[pair_3_y]==acids[pair_3_x]):
                observed+=1
        occurrences.append(observed)  


        fasta=0        
        for elem in background:
            if (elem[pair_1_y]==acids[pair_1_x]) and (elem[pair_2_y]==acids[pair_2_x]) and (elem[pair_3_y]==acids[pair_3_x]):
                fasta+=1

        expected=(fasta/len(background))*len(intervals)
        expactations.append(expected)

        if observed>10:
            chisq,p_value=chisquare([observed,len(intervals)-observed],f_exp=[expected,len(intervals)-expected])
            p_values.append(p_value)
        else:
            p_values.append('-')
            
    triple_motifs['Observed']=occurrences
    triple_motifs['Expected']=expactations
    triple_motifs['p-value']=p_values
    
    triple_motifs_selection=triple_motifs.where(triple_motifs['p-value']!='-').dropna()
    
    count=[]
    for i in range(len(triple_motifs_selection['p-value'].values)):
        if (triple_motifs_selection['p-value'].values[i])<0.05/len(triple_motifs['Location'].values):
            count.append(1)
        else:
            count.append(0)
    triple_motifs_selection['Count']=count
    
    triple_motifs_selection_p=triple_motifs_selection.where(triple_motifs_selection['Count']==1).dropna()
    triple_motifs_selection_p_copy=triple_motifs_selection_p.copy()
    del triple_motifs_selection_p_copy['Count']
    
    motifs=[]
    for i in range(len(triple_motifs_selection_p_copy['Location'].values)):
#        print(i)
        motif=dict()
        pair_1,pair_2,pair_3=triple_motifs_selection_p_copy['Location'].values[i]
        pair_1_x,pair_1_y=pair_1
        pair_2_x,pair_2_y=pair_2
        pair_3_x,pair_3_y=pair_3
        motif[pair_1_y]=acids[pair_1_x]
        motif[pair_2_y]=acids[pair_2_x]
        motif[pair_3_y]=acids[pair_3_x]
        motif[interval_length]=modification_site.lower()
        list_keys = list(motif.keys())
        list_keys.sort()
        word_motif=motif[list_keys[0]]+'.'*(list_keys[1]-list_keys[0]-1)+motif[list_keys[1]]+'.'*(
            list_keys[2]-list_keys[1]-1)+motif[list_keys[2]]+'.'*(list_keys[3]-list_keys[2]-1)+motif[list_keys[3]]
        motifs.append(word_motif)
        
    triple_motifs_selection_p_copy['Motifs']=motifs
    triple_motifs_selection_p_copy=triple_motifs_selection_p_copy.reset_index()    
    
    motifss=dict()
    counter=[]
    for i in range(len(triple_motifs_selection_p_copy['Location'].values)):
        if triple_motifs_selection_p_copy['Motifs'].values[i] not in motifss:
            motifss[triple_motifs_selection_p_copy['Motifs'].values[i]]=1
            counter.append(1)
        else:
            counter.append(0)
            
    triple_motifs_selection_p_copy['Counter']=counter
    triple_motifs_selection_p_copy_copy=triple_motifs_selection_p_copy.where(triple_motifs_selection_p_copy['Counter']==1).dropna()
    del triple_motifs_selection_p_copy_copy['Counter']
    del triple_motifs_selection_p_copy_copy['index']
    
    utils.saving_table(results_saving_dir,triple_motifs_selection_p_copy_copy,interval_length,'triple')
#    print('Triple motifs are ready!')
    logging.info(str(len(triple_motifs_selection_p_copy_copy['Motifs'].values))+u' triple (three acid length) motifs were identificated')
    
    return triple_motifs_selection_p_copy_copy



#делаем прогу для 4хбуквенного мотива
def quadruple_motifs(primary_motif,triple_motifs,acids,intervals,background,interval_length,results_saving_dir,modification_site):
    location=[]
    for i in range(len(triple_motifs['Location'].values)):
        pair_1,pair_2,pair_3=triple_motifs['Location'].values[i]
        for k in range(len(primary_motif['Location'].values)):
            pair_4=primary_motif['Location'].values[k]
            location.append((pair_1,pair_2,pair_3,pair_4))
            
    quadro_motifs=pd.DataFrame({'Location':location})
    
    count=[]
    for i in range(len(quadro_motifs['Location'].values)):
        pair_1,pair_2,pair_3,pair_4=quadro_motifs['Location'].values[i]
        pair_1_x,pair_1_y=pair_1
        pair_2_x,pair_2_y=pair_2
        pair_3_x,pair_3_y=pair_3
        pair_4_x,pair_4_y=pair_4
        if (pair_1_y!=pair_2_y) and (pair_1_y!=pair_3_y) and (pair_3_y!=pair_2_y) and (pair_4_y!=pair_2_y) and (pair_4_y!=pair_1_y) and (pair_4_y!=pair_3_y):
            count.append(1)
        else:
            count.append(0)
    quadro_motifs['Count']=count
    quadro_motifs=(quadro_motifs.where(quadro_motifs['Count']==1).dropna()).reset_index()
    del quadro_motifs['Count']
    del quadro_motifs['index']
    
    #для каждой тройки аминокислот должны рассчитать observed и expected
    occurrences=[]
    expactations=[]
    p_values=[]
    for i in range(len(quadro_motifs['Location'].values)):
#         print(i)
        pair_1,pair_2,pair_3,pair_4=quadro_motifs['Location'].values[i]
        pair_1_x,pair_1_y=pair_1
        pair_2_x,pair_2_y=pair_2
        pair_3_x,pair_3_y=pair_3
        pair_4_x,pair_4_y=pair_4
        observed=0
        for elem in intervals:
            if (elem[pair_1_y]==acids[pair_1_x]) and (elem[pair_2_y]==acids[pair_2_x]) and (elem[pair_3_y]==acids[pair_3_x]) and (elem[pair_4_y]==acids[pair_4_x]):
                observed+=1
        occurrences.append(observed)  


        fasta=0        
        for elem in background:
            if (elem[pair_1_y]==acids[pair_1_x]) and (elem[pair_2_y]==acids[pair_2_x]) and (elem[pair_3_y]==acids[pair_3_x]) and (elem[pair_4_y]==acids[pair_4_x]):
                fasta+=1

        expected=(fasta/len(background))*len(intervals)
        expactations.append(expected)

        if observed>10:
            chisq,p_value=chisquare([observed,len(intervals)-observed],f_exp=[expected,len(intervals)-expected])
            p_values.append(p_value)
        else:
            p_values.append('-')
            
    quadro_motifs['Observed']=occurrences
    quadro_motifs['Expected']=expactations
    quadro_motifs['p-value']=p_values
    
    quadro_motifs_selection=quadro_motifs.where(quadro_motifs['p-value']!='-').dropna()
    
    count=[]
    for i in range(len(quadro_motifs_selection['p-value'].values)):
        if (quadro_motifs_selection['p-value'].values[i])<0.05/len(quadro_motifs['Location'].values):
            count.append(1)
        else:
            count.append(0)
    quadro_motifs_selection['Count']=count
    
    quadro_motifs_selection_p=quadro_motifs_selection.where(quadro_motifs_selection['Count']==1).dropna()
    quadro_motifs_selection_p_copy=quadro_motifs_selection_p.copy()
    del quadro_motifs_selection_p_copy['Count']
    
    motifs=[]
    for i in range(len(quadro_motifs_selection_p_copy['Location'].values)):
#        print(i)
        motif=dict()
        pair_1,pair_2,pair_3,pair_4=quadro_motifs_selection_p_copy['Location'].values[i]
        pair_1_x,pair_1_y=pair_1
        pair_2_x,pair_2_y=pair_2
        pair_3_x,pair_3_y=pair_3
        pair_4_x,pair_4_y=pair_4
        motif[pair_1_y]=acids[pair_1_x]
        motif[pair_2_y]=acids[pair_2_x]
        motif[pair_3_y]=acids[pair_3_x]
        motif[pair_4_y]=acids[pair_4_x]
        motif[interval_length]=modification_site.lower()
        list_keys = list(motif.keys())
        list_keys.sort()
        word_motif=motif[list_keys[0]]+'.'*(list_keys[1]-list_keys[0]-1)+motif[list_keys[1]]+'.'*(
            list_keys[2]-list_keys[1]-1)+motif[list_keys[2]]+'.'*(list_keys[3]-list_keys[2]-1)+motif[list_keys[3]]+'.'*(
                list_keys[4]-list_keys[3]-1)+motif[list_keys[4]]
        motifs.append(word_motif)
        
    quadro_motifs_selection_p_copy['Motifs']=motifs
    quadro_motifs_selection_p_copy=quadro_motifs_selection_p_copy.reset_index()    
    
    motifss=dict()
    counter=[]
    for i in range(len(quadro_motifs_selection_p_copy['Location'].values)):
        if quadro_motifs_selection_p_copy['Motifs'].values[i] not in motifss:
            motifss[quadro_motifs_selection_p_copy['Motifs'].values[i]]=1
            counter.append(1)
        else:
            counter.append(0)
            
    quadro_motifs_selection_p_copy['Counter']=counter
    quadro_motifs_selection_p_copy_copy=quadro_motifs_selection_p_copy.where(quadro_motifs_selection_p_copy['Counter']==1).dropna()
    del quadro_motifs_selection_p_copy_copy['Counter']
    del quadro_motifs_selection_p_copy_copy['index']
    
    utils.saving_table(results_saving_dir,quadro_motifs_selection_p_copy_copy,interval_length,'quadruple')
    
    logging.info(str(len(quadro_motifs_selection_p_copy_copy['Motifs'].values))+u' quadruple (four acid length) motifs were identificated')
#    print('Quadro motifs are ready!')
    
    return quadro_motifs_selection_p_copy_copy


# In[117]:


#улучшенная версия 2
def motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,results_saving_dir):
    single=primary_motifs(acids,interval_length,chi2_selection,results_saving_dir,modification_site)
    double=double_motifs(single,acids,intervals,background,results_saving_dir,interval_length,modification_site)
    if double is not None:
        triple=triple_motifs(single,double,acids,intervals,background,interval_length,results_saving_dir,modification_site)
        if triple is not None:
            quadruple=quadruple_motifs(single,triple,acids,intervals,background,interval_length,results_saving_dir,modification_site)
        else:
            quadruple=None
    else:
        triple=None
    return single,double,triple,quadruple