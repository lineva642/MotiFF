# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:37:30 2020

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

def saving(working_dir,name_sample,interval_length,modification,modification_site,algorithm): 
    sample_saving_dir = os.path.join(working_dir,name_sample)
    results_saving_dir = os.path.join(sample_saving_dir,modification+'_'+
    modification_site+str(interval_length)+'_'+algorithm)
    if not (os.path.exists(sample_saving_dir)):
        os.mkdir(sample_saving_dir)
        if not (os.path.exists(results_saving_dir)):
            os.mkdir(results_saving_dir)
    else:
        if not (os.path.exists(results_saving_dir)):
            os.mkdir(results_saving_dir)
#    logging.basicConfig(level = logging.DEBUG,filename=os.path.join(results_saving_dir,'mylog.log'))        
#    logging.info(u'Directories for result saving are created')       
    return  sample_saving_dir,results_saving_dir


#поиск белков, соответствующих данным пептидам
def peptide_identification(path_FASTA,Peptides,sample_saving_dir):
#    print ('Identification of peptides')
    path_identificated_peptides=os.path.join(sample_saving_dir,'peptide_identification.csv')
    if not os.path.isfile(path_identificated_peptides):
        FASTA_dict=dict()
        for description, sequence in fasta.read(open(path_FASTA)):
            FASTA_dict[sequence]=description
        for elem in Peptides.index:
#            print(elem)
            #выделили элемент, будем искать его в FASTA
            for protein in FASTA_dict:
                if (protein.count(elem)!=0) and (Peptides['Protein'][elem] is not None):
                    Peptides['Protein'][elem]=None
                if (protein.count(elem)!=0) and (Peptides['Protein'][elem] is None):
                    Peptides['Protein'][elem]=protein            
        idPeptides=Peptides.dropna()
        
        idPeptides.to_csv(path_identificated_peptides,mode='w')
    else:
        idPeptides=pd.read_csv(path_identificated_peptides, sep=',')
        
#    print ('Peptides are identificated')
    logging.info(msg=str(len(idPeptides['Protein'].values))+' of '+str(len(Peptides['Protein'].values))
    +' modified peptides are identificated')
    logging.debug(msg=u'Peptides are identificated')                                
    return idPeptides


#нарезка интервалов
def interval_maker_experimental(idPeptides,interval_length,modification_site):
    intervals=[]
    #фиксируем пептид с модификациями
    #print('!',idPeptides)
    for elem in idPeptides.index:
        peptide,mod_peptide,protein = idPeptides['Peptide'][elem],idPeptides['Mod_Peptide'][elem],idPeptides['Protein'][elem]

        #число модификации в пептиде
        asterisk=0
        #бежим по выбранному пептиду
        for k in range(len(mod_peptide)):
            #проверяем, есть ли модификация
            if mod_peptide[k]=='*':

                pep_site=k-1-asterisk #позиция модификации в пептиде

                asterisk+=1


                #для каждой модификации нужно выделить интервал
                left=pep_site-interval_length
                right=pep_site+interval_length


                #возможно,этот интервал можно вырезать из нашего пептида, проверим эту гипотезу
                if ((left-interval_length)>=interval_length) and (right<=len(peptide)-1):
                    interval=peptide[left:right+1]
                    intervals.append(interval)
                #если так не работает, то нужно обратиться к белку
                else:
                    #добавим проверку на наличие символа '\n' в конце белка
                    if protein[-1]=='\n':
                        protein=protein[:-1]
                    a=protein.find(peptide)
                    protein_site=a+pep_site #позиция модификации в белке

                    left=protein_site-interval_length
                    right=protein_site+interval_length


                    #возможно,этот интервал можно вырезать из нашего белка, проверим эту гипотезу
                    if (left>=0) and (right<=len(protein)-1):
                        interval=protein[left:right+1]
                        intervals.append(interval)
                    
                    #рассматриваем случай, когда сайт прибит к N-концу белка
                    elif ((left<0) and (right<=len(protein)-1)):
                        interval=' '*(interval_length-protein_site)+protein[:right+1]
                        intervals.append(interval)
                        
                    #рассматриваем случай, когда сайт прибит к C-концу белка        
                    elif ((left>=0) and (right>len(protein)-1)):
                        interval=protein[left:]+' '*(right+1-len(protein))
                        intervals.append(interval)
                    
                    #рассматриваем случай, когда белок маленький 
                    else:
                        interval=' '*(interval_length-protein_site)+protein+' '*(right+1-len(protein))
                        intervals.append(interval)

    mod_intervals=[]

    for interval in intervals:
        if interval[interval_length]==modification_site:
            mod_intervals.append(interval)                   
    return mod_intervals


#формирование набора интервалов
def mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,sample_saving_dir):
#    print('Making intervals DB')
    idPeptides=peptide_identification(path_FASTA,Peptides,sample_saving_dir)
    mod_intervals=interval_maker_experimental(idPeptides,interval_length,modification_site)
#    print('Intervals DB is ready!')
    logging.debug(u'Intervals DB is ready')
    logging.info(u'Set of '+str(len(mod_intervals))+u' modified intervals s ready') 
    return mod_intervals

#правильная функция
def FASTA_DB_creation(path_FASTA):
    background_DB=dict()
    #берем белок из белковой БД человека
    for description, sequence in fasta.read(open(path_FASTA)):    
        #ищем его номер
        i=description.find('|')
        k=description.rfind('|')
        s=description[i+1:k]
        #если белка нет в background БД записываем 
        if s not in background_DB:
            if sequence.count('X')==0:
                background_DB[s]=sequence           
    return  background_DB    


#функция выделения интервалов для составления background
def Background_creator(elem,i,interval_length):
    #определим интервал поиска аминокислот        
    #рассматриваем случай, когда сайт прибит к N-концу белка
    if (i<=interval_length-1) and (len(elem)-i>interval_length):
        #количество кислот слева
        left=i
        #количество рассматриваемых кислот справа
        right=interval_length
        interval=' '*(interval_length-left)+elem[i-left:i+right+1]
    #рассматриваем случай, когда сайт прибит к C-концу белка    
    elif  (i>interval_length-1) and (len(elem)-i<=interval_length):
        left=interval_length
        right=len(elem)-i-1
        interval=elem[i-left:i+right+1]+' '*(interval_length-right)
    #рассматриваем случай, когда белок маленький    
    elif (i<=interval_length-1) and (len(elem)-i<=interval_length):
        left=i
        right=len(elem)-i-1
        interval=' '*(interval_length-left)+elem[i-left:i+right+1]+' '*(interval_length-right)
    #случай, когда можем рассмотреть интервал (-15,+15)    
    else:
        right=left=interval_length
        interval=elem[i-left:i+right+1]
    #выделяем наш интервал-последовательность из белка
    #interval=elem[number-left:number+right+1]
    return interval


#формирование набора интервалов из белоковой БД
def background_array_FASTA(path_FASTA,modification_site,interval_length):
    #хотим сделать background из идентифицированных белков
    background_DB=FASTA_DB_creation(path_FASTA)
#    print('background_array_FASTA',len(background_DB))
    background=[]

    for elem in background_DB:
        s=background_DB[elem]
        for i in range(len(s)):
            #ищем совпадающие с сайтом модификации
            if s[i]==modification_site:
                interval=Background_creator(s,i,interval_length)
                background.append(interval)               
    return background 


def background_maker(modification_site,interval_length,path_FASTA,modification):
#    print('Making background DB')
    background=background_array_FASTA(path_FASTA,modification_site,interval_length)
    logging.info(u'Set of ' + str(len(background))+ u' background intervals is created')
    logging.debug(u'Background DB is ready')
#    print('Background DB is ready')    
    return background   


#функции для валидации


#N_calculator= lambda intervals:len(intervals)


def n_calculator(acids,intervals,interval_length,results_saving_dir):
    #создаем матрицу для подсчета числа n и забиваем ее нулями
    n=[]
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        n.append(a)

    for interval in intervals:
        for k in range(len(interval)):
            #проверяем какая кислота из заданного набора стоит на зафиксированной позиции
            for i in range(len(acids)):
                #номер аминокислоты i в списке аминокислот совпадает с номером этой аминокислоты в матрице n
                acid=acids[i]
                #если совпадение есть,то добавляем единичку в матрицу n на соответствующую позицию
                if interval[k]==acid:
                    n[i][k]+=1
                else:
                    continue
    path=os.path.join(results_saving_dir,'occurrences.txt')
    saving=open(path,'w')
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            saving.write(str(n[i][k])+' ')
        saving.write('\n')
    saving.close()
    logging.debug(msg=u'Occurrence matrix was created')                    
    return n


def background_n_matrix(acids,interval_length,background,results_saving_dir):
    background_n=[]
    path=os.path.join(results_saving_dir,'background.txt')
    saving=open(path,'w')
    #забиваем матрицу нулями
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        background_n.append(a)
    #на данном этапе считаем встречаемость аминокислоты на j позиции
    for interval in background:
    #фиксируем позицию в интервале
        for k in range(len(interval)):
            #проверяем какая кислота из заданного набора стоит на зафиксированной позиции
            for i in range(len(acids)):
                #номер аминокислоты i в списке аминокислот совпадает с номером этой аминокислоты в матрице n
                acid=acids[i]
                #если совпадение есть,то добавляем единичку в матрицу n на соответствующую позицию
                if interval[k]==acid:
                    background_n[i][k]+=1
                else:
                    continue
              
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            if k!=interval_length:
                saving.write(str(background_n[i][k])+' ')
            else:
                saving.write(str(0)+' ')
        saving.write('\n')
    saving.close()
    logging.debug(msg=u'Background occurrence matrix was created')                     
    return background_n

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
    saving_table(results_saving_dir,primary_motif,interval_length,'primary')
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
    
    saving_table(results_saving_dir,second_motifs_selection_p_copy,interval_length,'double')
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
    
    saving_table(results_saving_dir,triple_motifs_selection_p_copy_copy,interval_length,'triple')
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
    
    saving_table(results_saving_dir,quadro_motifs_selection_p_copy_copy,interval_length,'quadruple')
    
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

##АЛГОРИТМ BINOMIAL

def p_value_bi(acids,modification_site,interval_length,background,background_n):
    background_n_1=[]
    for i in range(len(acids)):
        if acids[i]==modification_site:
            number=i
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        background_n_1.append(a)
    #составили матрицу вероятностей
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            if (k==interval_length) and (i==number):
                background_n_1[i][k]=1
            else:    
                background_n_1[i][k]=(background_n[i][k])/len(background)
    logging.debug(msg=u'p for each background occurrence matrix element was counted')            
    return background_n_1

def P_matrix_bi(acids,interval_length,n,N,background_n_1,results_saving_dir):
    P=[]
    path=os.path.join(results_saving_dir,'P.txt')
    saving=open(path,'w')
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        P.append(a)
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
#            print(i,k)
            #формула-binom.pmf(k, n, p, loc=0)
            result=0
            c=n[i][k]
            while c<=N:
                result=result+binom.pmf(c,N,background_n_1[i][k],loc=0)
                c+=1
            P[i][k]=result
            saving.write(str(result)+' ')
        saving.write('\n')
    saving.close()
    logging.debug(msg=u'Binomal probability for each occurrence matrix element was counted')    
    return P

def P_counter_bi(intervals,interval_length,modification_site,acids,background,results_saving_dir):
#    print('Probability is being counted')
    n=n_calculator(acids,intervals,interval_length,results_saving_dir)
    N=N_calculator(intervals)
    background_n=background_n_matrix(acids,interval_length,background,results_saving_dir)
    background_n_1=p_value_bi(acids,modification_site,interval_length,background,background_n)
    P=P_matrix_bi(acids,interval_length,n,N,background_n_1,results_saving_dir)
#    print('Probability is being counted')
    logging.info(u'Binomial probability for each amino acid in matrix was counted')
    return P

#рассчитаем occurrences
def occurrences_counter_bi(intervals,interval_length,acids,modification_site,results_saving_dir):
    n=n_calculator(acids,intervals,interval_length,results_saving_dir)
#    N=N_calculator(intervals)
    path=os.path.join(results_saving_dir,'occurrences.txt')
    saving=open(path,'w')
    occurrences=[]
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        occurrences.append(a)
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            if k!=interval_length:
                occurrences[i][k]=n[i][k]
            saving.write(str(occurrences[i][k])+' ')
        saving.write('\n')
    saving.close()             
    logging.debug(u'Occurrences for each amino acid in matrix was counted')            
    return occurrences 

#отбор по вероятности
def P_validation_bi(acids,interval_length,P):
    P1=[]
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        P1.append(a)
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            if P[i][k]<0.05/(21*13):
                continue
            else:
                P1[i][k]+=1
    logging.debug(u'Occurrences for each amino acid in matrix was counted')            
    return P1

#отбор по встречаемости более 10
def occurrences_validation_bi(acids,interval_length,occurrences):
    occurrences_1=[]
    for i in range(len(acids)):
        a=[0]*(interval_length*2+1)
        occurrences_1.append(a)
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            if occurrences[i][k]>=10:
                occurrences_1[i][k]=occurrences[i][k]
    return occurrences_1

#соединим валидацию по occurances и p-values
#нужно совместить occurances_1 & P1
def final_validation_bi(acids,interval_length,occurrences,P):
    occurrences_1=occurrences_validation_bi(acids,interval_length,occurrences)
    P1=P_validation_bi(acids,interval_length,P)
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            if occurrences_1[i][k]>=10:
                continue
            else:
                P1[i][k]=1
    (u'Binomial probabilities were counted and selection with p=0.05/Bonferroni correction and occurrences>10 was performed')            
    return P1 


def single_motifs_creator_bi(acids,final_P,P,background,intervals,interval_length,modification_site):
    indexes=[]
    location=[]
    probability=[]
    occurrences=[]
    backgrounds=[]
    for i in range(len(final_P[:])):
        for j in range(len(final_P[0])):
            if final_P[i][j]==0 and acids[i]!=' ':
                #расчет индексов
                if j<interval_length:
                    s = [acids[i], '.'*(interval_length-j-1), modification_site.lower()]
                    index=''.join(s)
                    #index=acids[i]+'.'*(interval_length-j-1)+modification_site.lower()
                else:
                    s = [modification_site.lower(), '.'*(j-interval_length-1), acids[i]]
                    index=''.join(s)                    
                    #index=modification_site.lower()+'.'*(j-interval_length-1)+acids[i]
                x=sum(1 for interval in background if ((interval[j]==acids[i])))
                u=sum(1 for interval in intervals if ((interval[j]==acids[i])))
                indexes.append(index)    
                location.append((i,j))
                probability.append(P[i][j])
                occurrences.append(u)
                backgrounds.append(x)
    #создадим Series
    single_motifs=pd.DataFrame({'Location': location, 'Probability': probability , 'Occurrences' : occurrences, 'Backgrounds': backgrounds}, index=indexes)
    logging.info(str(len(single_motifs['Occurrences'].values))+u' primary (one acid length) motifs were identificated')    
    return single_motifs 


def double_motifs_creator_bi(acids,single_motifs_creator,background,intervals,P,interval_length,modification_site,results_saving_dir):
    #для каждого претендента нужно построить двойной комплексный мотив и оценить вероятность такого мотива
    indexes=[]
    location=[]
    probability=[]
    occurrences=[]
    backgrounds=[]
    P_double=[]
    P_B_double=[]

        
    #закрепили одну аминокислоту
    for i in range(len(single_motifs_creator)):
        #выбрали вторую аминокислоту
        #print(i)
        for j in range(len(single_motifs_creator)):
            if j!=i:
                #print(j)
                #выбрали наш сложный мотив как две АА
                motif_1_x,motif_1_y=single_motifs_creator['Location'].values[i]
                motif_2_x,motif_2_y=single_motifs_creator['Location'].values[j]

                #рассматриваем кислоты, которые стоят в интервале на разных позициях
                if motif_1_y!=motif_2_y:
                    
                    #считаем количество интервалов с данным комплексным мотивом в background
                    x=sum(1 for interval in background if ((interval[motif_1_y]==acids[motif_1_x]) and (interval[motif_2_y]==acids[motif_2_x])))

                    #считаем p-value этого комплексного мотива
                    p=x/len(background)

                    #теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
                    n=sum(1 for interval in intervals if ((interval[motif_1_y]==acids[motif_1_x]) and (interval[motif_2_y]==acids[motif_2_x])))

                    #считаем Р вероятность по биномиальной формуле для комплексного мотива
                    P_AB=0
                    c=n
                    while c<=len(intervals):
                        P_AB=P_AB+binom.pmf(n,len(intervals),p,loc=0)
                        c+=1
                    #print('P_AB:',P_AB)

                    #нашли вероятность P(AB), теперь найдем условную вероятность
                    P_B=P[motif_1_x][motif_1_y]
                    #print('P_B:',P_B,motif_1_x,motif_1_y)
                    P_A_B=P_AB/P_B
                    #print('P_A_B:',P_A_B)

                    
                    #создадим dataframe
                    if (motif_1_y<motif_2_y) and (motif_1_y<interval_length) and (motif_2_y>interval_length):
                        s = [acids[motif_1_x], '.'*(interval_length-motif_1_y-1), modification_site.lower(),'.'*(motif_2_y-interval_length-1),acids[motif_2_x]]
                        index=''.join(s)
                        #s=motif_1_acid+'.'*(interval_length-motif_1_y-1)+modification_site.lower()+'.'*(motif_2_y-interval_length-1)+motif_2_acid
                        indexes.append(index)
                    elif (motif_1_y<motif_2_y) and (motif_1_y<interval_length) and (motif_2_y<interval_length):
                        s = [acids[motif_1_x], '.'*(motif_2_y-motif_1_y-1),acids[motif_2_x],'.'*(interval_length-motif_2_y-1),modification_site.lower()]
                        index=''.join(s)
                        #s=motif_1_acid+'.'*(motif_2_y-motif_1_y-1)+motif_2_acid+'.'*(interval_length-motif_2_y-1)+modification_site.lower()
                        indexes.append(index)
                    elif (motif_1_y<motif_2_y) and (motif_1_y>interval_length):
                        s = [modification_site.lower(), '.'*(motif_1_y-interval_length-1),acids[motif_1_x],'.'*(motif_2_y-motif_1_y-1),acids[motif_2_x]]
                        index=''.join(s)
                        #s=modification_site.lower()+'.'*(motif_1_y-interval_length-1)+motif_1_acid+'.'*(motif_2_y-motif_1_y-1)+motif_2_acid
                        indexes.append(index)
                    elif (motif_2_y<motif_1_y) and (motif_2_y<interval_length) and (motif_1_y>interval_length):
                        s = [acids[motif_2_x], '.'*(interval_length-motif_2_y-1),modification_site.lower(),'.'*(motif_1_y-interval_length-1),acids[motif_1_x]]
                        index=''.join(s)
                        #s=motif_2_acid+'.'*(interval_length-motif_2_y-1)+modification_site.lower()+'.'*(motif_1_y-interval_length-1)+motif_1_acid
                        indexes.append(index)
                    elif (motif_2_y<motif_1_y) and (motif_2_y<interval_length) and (motif_1_y<interval_length):
                        s = [acids[motif_2_x], '.'*(motif_1_y-motif_2_y-1),acids[motif_1_x],'.'*(interval_length-motif_1_y-1),modification_site.lower()]
                        index=''.join(s)                        
                        #s=motif_2_acid+'.'*(motif_1_y-motif_2_y-1)+motif_1_acid+'.'*(interval_length-motif_1_y-1)+modification_site.lower()
                        indexes.append(index)
                    elif (motif_2_y<motif_1_y) and (motif_2_y>interval_length):
                        s = [modification_site.lower(), '.'*(motif_2_y-interval_length-1),acids[motif_2_x],'.'*(motif_1_y-motif_2_y-1),acids[motif_1_x]]
                        index=''.join(s)                            
                        #s=modification_site.lower()+'.'*(motif_2_y-interval_length-1)+motif_2_acid+'.'*(motif_1_y-motif_2_y-1)+motif_1_acid
                        indexes.append(index)

                    if motif_1_y<motif_2_y:
                        location.append(((motif_1_x,motif_1_y),(motif_2_x,motif_2_y)))
                    else:
                        location.append(((motif_2_x,motif_2_y),(motif_1_x,motif_1_y)))
                    #location.append((motif_2_x,motif_2_y))

                    probability.append(P_A_B)

                    occurrences.append(n)
                    
                    backgrounds.append(x)

                    P_double.append(P_AB)
                    
                    P_B_double.append(P_B)
                    
                    
    motifs=pd.DataFrame({'Location': location, 'Probability': probability , 'Occurrences' : occurrences, 'Backgrounds': backgrounds, 'P_AB' : P_double, 'P_B' : P_B_double}, index=indexes)
    sorted_motifs=motifs.copy()
    sorted_motifs.sort_values('Probability',ascending=True,inplace=True)
    double_motifs=sorted_motifs
    
    count=[]
    for i in range(len(double_motifs['Probability'].values)):
        if ((double_motifs['Probability'].values[i])<0.05/len(double_motifs['Probability'].values) and (double_motifs['Occurrences'].values[i])>10):
            count.append(1)
        else:
            count.append(0)
    double_motifs['Count']=count
    
    double_motifs_selection=double_motifs.where(double_motifs['Count']==1).dropna()
    double_motifs_selection_copy=double_motifs_selection.copy()
    del double_motifs_selection_copy['Count']
    
    
#    print('double_motifs',double_motifs_selection_copy)
    
    if probability==[]:
        double_motifs_selection_copy=None
    
#     if double_motifs_selection_copy is not None:
#         volcano_plot_for_motifs(double_motifs_selection_copy,path_results,name='double_motifs')    
        
    logging.info(str(len(double_motifs_selection_copy['Occurrences'].values))+u' double (two acid length) motifs were identificated')
    return double_motifs_selection_copy

def triple_motifs_creator_bi(acids,double_motifs,single_motifs,background,intervals,modification_site,interval_length,results_saving_dir):

    indexes=[]
    location=[]
    probability=[]
    occurrences=[]
    backgrounds=[]
    P_triple=[]
    P_BC_triple=[]
    #зафиксировали двубуквенную часть
    for k in range(len(double_motifs['Location'].values)):
        elem=double_motifs['Location'].values[k]
        double_motif_1,double_motif_2=elem
        double_motif_1_x,double_motif_1_y=double_motif_1
        double_motif_2_x,double_motif_2_y=double_motif_2
        
        #теперь хотим добавлять однобуквенную 
        #double_motif_AA=double_motifs['Location'].index[k]
        for i in range(len(single_motifs['Location'].values)):
            third_x,third_y=single_motifs['Location'].values[i]
            #third_AA=acids[third_x]
            #порверяем, что выбранная кислота не совпадает с предыдущими и нет накладок
            if (third_y!=double_motif_1_y) and (third_y!=double_motif_2_y):
                
                    #считаем количество интервалов с тройным мотивом в background
                    x=sum(1 for interval in background if ((interval[double_motif_1_y]==acids[double_motif_1_x]) 
                                                           and (interval[double_motif_2_y]==acids[double_motif_2_x]) 
                                                           and (interval[third_y]==acids[third_x])))

                    #считаем p-value этого комплексного мотива
                    p=x/len(background)

                    #теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
                    n=sum(1 for interval in intervals if ((interval[double_motif_1_y]==acids[double_motif_1_x]) 
                                                           and (interval[double_motif_2_y]==acids[double_motif_2_x]) 
                                                           and (interval[third_y]==acids[third_x])))

                    #считаем Р вероятность по биномиальной формуле для тройного мотива
                    P_ABC=0
                    c=n
                    while c<=len(intervals):
                        P_ABC=P_ABC+binom.pmf(n,len(intervals),p,loc=0)
                        c+=1
                    #print('P_AB:',P_AB)

                    #нашли вероятность P(ABC), теперь найдем условную вероятность
                    P_BC=double_motifs['P_AB'].values[k]      
                    #print('P_B:',P_B,motif_1_x,motif_1_y)
                    P_A_BC=P_ABC/P_BC
                    #print('P_A_B:',P_A_B)
                
                    #теперь разбираемся с названием мотива
                    names=dict()
                    names[double_motif_1_y]=acids[double_motif_1_x]
#                    print('names dict',double_motif_1_y,acids[double_motif_1_x])
                    names[double_motif_2_y]=acids[double_motif_2_x]
#                    print('names dict',double_motif_2_y,acids[double_motif_2_x])
                    names[third_y]=acids[third_x]
#                    print('names',third_y,acids[third_x])
                    names[interval_length]=modification_site.lower()
                    list_names=list(names.keys())
#                    print('list names 0',list_names)
                    list_names.sort()
#                    print('list names',list_names)
                    
                    name=[names[list_names[0]],'.'*(list_names[1]-list_names[0]-1),names[list_names[1]],
                          '.'*(list_names[2]-list_names[1]-1),names[list_names[2]],'.'*(list_names[3]-list_names[2]-1),
                          names[list_names[3]]]
                    motif=''.join(name)
                          
                    loc=((acids.index(names[list_names[0]].upper()),list_names[0]), (acids.index(names[list_names[1]].upper()),list_names[1]),
                          (acids.index(names[list_names[2]].upper()),list_names[2]),(acids.index(names[list_names[3]].upper()),list_names[3]))
                          
                    probability.append(P_A_BC)

                    occurrences.append(n)
                    
                    backgrounds.append(x)

                    P_triple.append(P_ABC)
                    
                    P_BC_triple.append(P_BC)
                          
                    location.append(loc)
                          
                    indexes.append(motif)      
                    
                    
    motifs=pd.DataFrame({'Location': location, 'Probability': probability , 'Occurrences' : occurrences, 'Backgrounds':backgrounds, 'P_ABC' : P_triple, 'P_BC' : P_BC_triple}, index=indexes)
    sorted_motifs=motifs.copy()
    sorted_motifs.sort_values('Probability',ascending=True,inplace=True)
    triple_motifs=sorted_motifs
    
    count=[]
    for i in range(len(triple_motifs['Probability'].values)):
        if ((triple_motifs['Probability'].values[i])<0.05/len(triple_motifs['Probability'].values) and (triple_motifs['Occurrences'].values[i])>10):
            count.append(1)
        else:
            count.append(0)
    triple_motifs['Count']=count
    
    triple_motifs_selection=triple_motifs.where(triple_motifs['Count']==1).dropna()
    triple_motifs_selection_copy=triple_motifs_selection.copy()
    del triple_motifs_selection_copy['Count']
       
    if probability==[]:
        triple_motifs_selection_copy=None
    
#     if triple_motifs_selection_copy is not None:
# #         volcano_plot_for_motifs(triple_motifs_selection_copy,path_results,name='triple_motifs')    
    logging.info(str(len(triple_motifs_selection_copy['Occurrences'].values))+u' triple (three acid length) motifs were identificated')    
    return triple_motifs_selection_copy

def quadruple_motifs_creator_bi(acids,triple_motifs,single_motifs,background,intervals,modification_site,interval_length,results_saving_dir):

    indexes=[]
    location=[]
    probability=[]
    occurrences=[]
    backgrounds=[]
    P_triple=[]
    P_BCD_triple=[]
    #зафиксировали трехбуквенную часть
    for k in range(len(triple_motifs['Location'].values)):
        elem=triple_motifs['Location'].values[k]
        triple_motif_1,triple_motif_2,triple_motif_3,triple_motif_4=elem
        triple_motif_1_x,triple_motif_1_y=triple_motif_1
        triple_motif_2_x,triple_motif_2_y=triple_motif_2
        triple_motif_3_x,triple_motif_3_y=triple_motif_3
        triple_motif_4_x,triple_motif_4_y=triple_motif_4
        
        #теперь хотим добавлять однобуквенную 
        #double_motif_AA=double_motifs['Location'].index[k]
        for i in range(len(single_motifs['Location'].values)):
            fourth_x,fourth_y=single_motifs['Location'].values[i]
            #third_AA=acids[third_x]
            #порверяем, что выбранная кислота не совпадает с предыдущими и нет накладок
            if (fourth_y!=triple_motif_1_y) and (fourth_y!=triple_motif_2_y) and (fourth_y!=triple_motif_3_y) and (fourth_y!=triple_motif_4_y):
                
                    #считаем количество интервалов с тройным мотивом в background
                    x=sum(1 for interval in background if ((interval[triple_motif_1_y]==acids[triple_motif_1_x]) 
                                                           and (interval[triple_motif_2_y]==acids[triple_motif_2_x]) 
                                                           and (interval[fourth_y]==acids[fourth_x]))
                                                           and (interval[triple_motif_3_y]==acids[triple_motif_3_x])
                                                           and (interval[triple_motif_4_y]==acids[triple_motif_4_x])) 
                                                        

                    #считаем p-value этого комплексного мотива
                    p=x/len(background)

                    #теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
                    n=sum(1 for interval in intervals if ((interval[triple_motif_1_y]==acids[triple_motif_1_x]) 
                                                           and (interval[triple_motif_2_y]==acids[triple_motif_2_x]) 
                                                           and (interval[fourth_y]==acids[fourth_x]))
                                                           and (interval[triple_motif_3_y]==acids[triple_motif_3_x])
                                                           and (interval[triple_motif_4_y]==acids[triple_motif_4_x]))
                    #считаем Р вероятность по биномиальной формуле для тройного мотива
                    P_ABCD=0
                    c=n
                    while c<=len(intervals):
                        P_ABCD=P_ABCD+binom.pmf(n,len(intervals),p,loc=0)
                        c+=1
                    #print('P_AB:',P_AB)

                    #нашли вероятность P(ABC), теперь найдем условную вероятность
                    P_BCD=triple_motifs['P_ABC'].values[k]      
                    #print('P_B:',P_B,motif_1_x,motif_1_y)
                    P_A_BCD=P_ABCD/P_BCD
                    #print('P_A_B:',P_A_B)
                
                    #теперь разбираемся с названием мотива
                    names=dict()
                    y=[triple_motif_1_y,triple_motif_2_y,triple_motif_3_y,triple_motif_4_y,fourth_y]
                    x=[triple_motif_1_x,triple_motif_2_x,triple_motif_3_x,triple_motif_4_x,fourth_x]
#                    print(len(y))
                    for i in range(len(y)):
                        igrek=y[i]
                        ixes=x[i]
                        if igrek!=interval_length:
                            names[igrek]=acids[ixes]
                        else:
                            names[igrek]=modification_site.lower()
                    
                    list_names=list(names.keys())
#                    print('list names 0',list_names)
                    list_names.sort()
#                    print('list names',list_names)
                    
                    name=[names[list_names[0]],'.'*(list_names[1]-list_names[0]-1),names[list_names[1]],
                          '.'*(list_names[2]-list_names[1]-1),names[list_names[2]],'.'*(list_names[3]-list_names[2]-1),
                          names[list_names[3]],'.'*(list_names[4]-list_names[3]-1),names[list_names[4]]]
                    motif=''.join(name)
                          
                    loc=((acids.index(names[list_names[0]].upper()),list_names[0]), (acids.index(names[list_names[1]].upper()),list_names[1]),
                          (acids.index(names[list_names[2]].upper()),list_names[2]),(acids.index(names[list_names[3]].upper()),list_names[3]),
                          (acids.index(names[list_names[4]].upper()),list_names[4]))
                          
                    probability.append(P_A_BCD)

                    occurrences.append(n)
                    
                    backgrounds.append(x)

                    P_triple.append(P_ABCD)
                    
                    P_BCD_triple.append(P_BCD)
                          
                    location.append(loc)
                          
                    indexes.append(motif)      
                    
                    
    motifs=pd.DataFrame({'Location': location, 'Probability': probability , 'Occurrences' : occurrences, 'Backgrounds':backgrounds, 'P_ABCD' : P_triple, 'P_BCD' : P_BCD_triple}, index=indexes)
    sorted_motifs=motifs.copy()
    sorted_motifs.sort_values('Probability',ascending=True,inplace=True)
    fourth_motifs=sorted_motifs
    
    count=[]
    for i in range(len(fourth_motifs['Probability'].values)):
        if ((fourth_motifs['Probability'].values[i])<0.05/len(fourth_motifs['Probability'].values) and (fourth_motifs['Occurrences'].values[i])>10):
            count.append(1)
        else:
            count.append(0)
    fourth_motifs['Count']=count
    
    fourth_motifs_selection=fourth_motifs.where(fourth_motifs['Count']==1).dropna()
    fourth_motifs_selection_copy=fourth_motifs_selection.copy()
    del fourth_motifs_selection_copy['Count']
       
    if probability==[]:
        fourth_motifs_selection_copy=None
    
#     if fourth_motifs_selection_copy is not None:
#         volcano_plot_for_motifs(fourth_motifs,path_results,name='triple_motifs')    
    logging.info(str(len(fourth_motifs_selection_copy['Occurrences'].values))+u' quadruple (four acid length) motifs were identificated')    
    return fourth_motifs_selection_copy

def motifs_bi(acids,P_final,interval_length,modification_site,background,intervals,P,results_saving_dir):
    #primary=primary_motifs_creator(acids,P_final,interval_length,modification_site)
    single=single_motifs_creator_bi(acids,P_final,P,background,intervals,interval_length,modification_site)
    saving_table(results_saving_dir,single,interval_length,'single')
    double=double_motifs_creator_bi(acids,single,background,intervals,P,interval_length,modification_site,results_saving_dir)
    if double is not None:
        saving_table(results_saving_dir,double,interval_length,'double')
        triple=triple_motifs_creator_bi(acids,double,single,background,intervals,modification_site,interval_length,results_saving_dir)
    if triple is not None:
        saving_table(results_saving_dir,triple,interval_length,'triple')
        quadruple=quadruple_motifs_creator_bi(acids,triple,single,background,intervals,modification_site,interval_length,results_saving_dir)
    if quadruple is not None:
        saving_table(results_saving_dir,quadruple,interval_length,'quadruple')
    return single,double,triple,quadruple

def output_experimental_bi(working_dir,name_sample,Peptides,interval_length,modification_site,modification,path_FASTA,algorithm):
    sample_saving_dir,results_saving_dir=saving(working_dir,name_sample,interval_length,modification,modification_site,algorithm)
    acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
    intervals=mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,sample_saving_dir)
#    print('Длина интервалов в exp',len(intervals))
    background=background_maker(modification_site,interval_length,path_FASTA,modification)
#    print('Длина интервалов в back',len(background))
    P=P_counter_bi(intervals,interval_length,modification_site,acids,background,results_saving_dir)
    occurrences=occurrences_counter_bi(intervals,interval_length,acids,modification_site,results_saving_dir)
    P_final=final_validation_bi(acids,interval_length,occurrences,P)
    single,double,triple,quadruple=motifs_bi(acids,P_final,interval_length,modification_site,background,intervals,P,results_saving_dir)
#    print(single,double,triple,quadruple)
#    logging.info(msg='Program was finished successfully')
    return P,occurrences,intervals,background

def saving_table(results_saving_dir,result,interval_length,name):
    path=os.path.join(results_saving_dir,'table'+str(interval_length)+'_'+name+'.csv')
    result.to_csv(path)   


# Для алгоритма chi2


def output_experimental(working_dir,name_sample,Peptides,interval_length,modification_site,modification,path_FASTA):
    #path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
    #path_identification='/home/vikalin/Article/identification.txt'
    sample_saving_dir,results_saving_dir=saving(working_dir,name_sample,interval_length,modification,modification_site,'chi2')
#    path_sample,path_modification,path_results=saving(name_sample,interval_length,modification,modification_site,'chi2')
    acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
    intervals=mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,sample_saving_dir)
#    print('Длина интервалов в exp',len(intervals))
    background=background_maker(modification_site,interval_length,path_FASTA,modification)
#    print('Длина интервалов в back',len(background))
    occurrences=n_calculator(acids,intervals,interval_length,results_saving_dir)
    background_n=background_n_matrix(acids,interval_length,background,results_saving_dir)
    expected_FASTA,chi2_results,chi2_selection=p_value(background_n,occurrences,interval_length,modification_site,acids,results_saving_dir)
    single,double,triple,quadruple=motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,results_saving_dir)
#    print(single,double,triple,quadruple)
#    logging.info(msg='Program was finished successfully')

    
    return chi2_results,chi2_selection,intervals,background

def output(algorithm,working_dir,name_sample,Peptides,interval_length,modification_site,modification,path_FASTA):
    if algorithm=="binom":
#        print('I AM HERE')
        a,b,c,d=output_experimental_bi(working_dir,name_sample,Peptides,interval_length,
                                       modification_site,modification,path_FASTA,algorithm)
        #a,b,c,d=P,occurrences,intervals,background
    else:
        a,b,c,d=output_experimental(working_dir,name_sample,Peptides,interval_length,
                                    modification_site,modification,path_FASTA)
        #a,b,c,d=chi2_results,chi2_selection,intervals,background
    logging.info(msg='Program was finished successfully')                                
    return  a,b,c,d  
