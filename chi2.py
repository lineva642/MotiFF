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


#def p_value_selection(chi2_results,occurrences_FASTA,interval_length):
#    chi2_selection=[]
#    for i in range(0,21):
#        a=[0]*(2*interval_length+1)
#        chi2_selection.append(a)
#    for i in range(0,21):
#        for k in range(0,2*interval_length+1):
#            if (chi2_results[i][k]<0.05/(21*(2*interval_length+1))) and (k!=interval_length) and (occurrences_FASTA[i][k]>10):
#                chi2_selection[i][k]=1
#    logging.debug(u'Selection of matrix elements was successful')            
#    return chi2_selection            

def p_value_selection(chi2_results,occurrences_FASTA,interval_length):
    chi2_selection=[]
#    for i in range(0,21):
#        a=[0]*(2*interval_length+1)
#        chi2_selection.append(a)
    for i in range(0,21):
        for k in range(0,2*interval_length+1):
            if (chi2_results[i][k]<0.05/(21*(2*interval_length+1))) and (k!=interval_length) and (occurrences_FASTA[i][k]>10):
                chi2_selection.append((i,k))
    logging.debug(u'Selection of matrix elements was successful')          
    return chi2_selection            

# In[108]:


#напишем одну функцию для подсчета p-value
def p_value(background_n,occurrences,interval_length,modification_site,acids,results_saving_dir):
    expected_FASTA,all_exp=expected_destribution(background_n,occurrences,interval_length)
    chi2_results=chi2_result(occurrences,expected_FASTA,all_exp,interval_length,results_saving_dir)
    chi2_selection=p_value_selection(chi2_results,occurrences,interval_length)
#    heatmap_visualization(chi2_selection,acids,interval_length,modification_site,'Отбор по p-value',results_saving_dir,name='p_value_selection')
    logging.info(u'Chi2 values were counted and selection with p=0.05/Bonferroni and occurrences>10 correction was performed')
    return expected_FASTA,chi2_results,chi2_selection

#функция для подсчета комлексных мотивов
def primary_motifs(acids,interval_length,chi2_selection,results_saving_dir,modification_site):    
    primary_motifs_number=[]
    primary_motifs_letter=[]
    for elem in chi2_selection:
        i,k=elem
        motif=np.zeros(interval_length*2+1)
        motif[k]=i
        primary_motifs_number.append(motif)                    
        if k<interval_length:
            motif=acids[i]+'.'*(interval_length-k-1)+modification_site.lower()
        else:
            motif=modification_site.lower()+'.'*(k-interval_length-1)+acids[i]
        primary_motifs_letter.append(motif)
    vector=np.array(primary_motifs_number)    
    table=pd.DataFrame({'Number motif':primary_motifs_number,'Letter motif':primary_motifs_letter})      
    utils.saving_table(results_saving_dir,table,interval_length,'primary')
    logging.info(str(len(table['Number motif'].values))+u' primary (one acid length) motifs were identificated')
    return vector,table


# In[114]:
def counter_doubles(dataset,i_1,i_2,k_1,k_2,acids):
    count=0
    AA_1,AA_2=acids[i_1],acids[i_2]
    for elem in dataset:
        if (elem[k_1]==AA_1) and (elem[k_2]==AA_2):
            count+=1
    return count

def double_motifs(vector,acids,intervals,background,results_saving_dir,interval_length,modification_site):
    #выполняем тензорное умножение
    b = np.tensordot(vector, vector.T, axes=0)
    l=len(vector)    
    matrix=np.zeros((l,l,interval_length*2+1))
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
                observed=counter_doubles(intervals,i_1,i_2,k_1,k_2,acids)
                fasta=counter_doubles(background,i_1,i_2,k_1,k_2,acids)
                expected=(fasta/back_len)*int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])
                position={k_1:acids[i_1],k_2:acids[i_2],interval_length:modification_site.lower()}
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

def counter_triples(dataset,i_1,i_2,i_3,k_1,k_2,k_3,acids):
    count=0
    AA_1,AA_2,AA_3=acids[i_1],acids[i_2],acids[i_3]
    for elem in dataset:
        if (elem[k_1]==AA_1) and (elem[k_2]==AA_2) and (elem[k_3]==AA_3):
            count+=1
    return count

def triple_motifs(primary_motifs,double_motifs,acids,intervals,background,interval_length,results_saving_dir,modification_site):
    #составили всевозможные пары тройных мотивов
#    table=((double_motifs.where(double_motifs['Observed']>=10)).dropna()).reset_index()
    b=len(double_motifs['Observed'].values)
    table=(double_motifs.loc[double_motifs['p-value']<0.05/b][double_motifs['Observed']>=10]).reset_index()
    del table['index']

    double_vector=np.zeros((len(table['Number motif'].values),interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        double_vector[i,:]=table['Number motif'].values[i]
    print(double_vector[0])

    vector=np.zeros((len(primary_motifs['Number motif'].values),interval_length*2+1)) 
    for i,ik in enumerate(primary_motifs['Number motif'].values):
        vector[i,:]=primary_motifs['Number motif'].values[i]
    print(vector.shape)
    
    b=np.tensordot(double_vector,vector.T,axes=0)
    matrix=np.zeros((len(vector),len(double_vector),interval_length*2+1))
    back_len=len(background)
    int_len=len(intervals)
    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0]):
            matrix[j,i,:]=np.diag(b[i,:,:,j])
            
#    print(matrix.shape)
#    print(matrix[0][0],vector[0],double_vector[0])
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})     
    for i in range(len(vector)):
        for j in range(len(double_vector)):
            elem=matrix[i,j]
            print('elem',elem)
            if (elem.any())==False:
                motif=vector[i]+double_vector[j]
                print('motif',motif)
                #положение аминокислот в интервале
                
                k_1,k_2,k_3=np.nonzero(motif)[0][0],np.nonzero(motif)[0][1],np.nonzero(motif)[0][2]
                #номера аминокислот
                i_1,i_2,i_3=int(motif[k_1]),int(motif[k_2]),int(motif[k_3])
                print('!',k_1,k_2,k_3,i_1,i_2,i_3)
                observed=counter_triples(intervals,i_1,i_2,i_3,k_1,k_2,k_3,acids)
                fasta=counter_triples(background,i_1,i_2,i_3,k_1,k_2,k_3,acids)
                expected=(fasta/back_len)*int_len
                chisq,p_value=chisquare([observed,int_len-observed],f_exp=[expected,int_len-expected])
                position={k_1:acids[i_1],k_2:acids[i_2],k_3:acids[i_3],interval_length:modification_site.lower()}
                keys=list(position.keys())
                keys.sort()
                motif_l=''.join([position[keys[0]],'.'*(keys[1]-keys[0]-1),position[keys[1]],'.'*(keys[2]-keys[1]-1),position[keys[2]],'.'*(keys[3]-keys[2]-1),position[keys[3]]])
                result=result.append({'Letter motif':motif_l,'Number motif':motif,
                                'Observed':observed,'Expected':expected,
                                                    'p-value':p_value},
                                                    ignore_index=True)               
            
            
            
    return result


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
    vector,single=primary_motifs(acids,interval_length,chi2_selection,results_saving_dir,modification_site)
    double=double_motifs(vector,acids,intervals,background,results_saving_dir,interval_length,modification_site)
    if double is not None:
        triple=triple_motifs(single,double,acids,intervals,background,interval_length,results_saving_dir,modification_site)
#        if triple is not None:
#            quadruple=quadruple_motifs(single,triple,acids,intervals,background,interval_length,results_saving_dir,modification_site)
#        else:
#            quadruple=None
#    else:
#        triple=None
#    return single,double,triple,quadruple
    return single,double,triple