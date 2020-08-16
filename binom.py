# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:55:43 2020

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
    n=utils.n_calculator(acids,intervals,interval_length,results_saving_dir)
    N_calculator= lambda intervals:len(intervals)
    N=N_calculator(intervals)
    background_n=utils.background_n_matrix(acids,interval_length,background,results_saving_dir)
    background_n_1=p_value_bi(acids,modification_site,interval_length,background,background_n)
    P=P_matrix_bi(acids,interval_length,n,N,background_n_1,results_saving_dir)
#    print('Probability is being counted')
    logging.info(u'Binomial probability for each amino acid in matrix was counted')
    return P

#рассчитаем occurrences
def occurrences_counter_bi(intervals,interval_length,acids,modification_site,results_saving_dir):
    n=utils.n_calculator(acids,intervals,interval_length,results_saving_dir)
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
    utils.saving_table(results_saving_dir,single,interval_length,'single')
    double=double_motifs_creator_bi(acids,single,background,intervals,P,interval_length,modification_site,results_saving_dir)
    if double is not None:
        utils.saving_table(results_saving_dir,double,interval_length,'double')
        triple=triple_motifs_creator_bi(acids,double,single,background,intervals,modification_site,interval_length,results_saving_dir)
    if triple is not None:
        utils.saving_table(results_saving_dir,triple,interval_length,'triple')
        quadruple=quadruple_motifs_creator_bi(acids,triple,single,background,intervals,modification_site,interval_length,results_saving_dir)
    if quadruple is not None:
        utils.saving_table(results_saving_dir,quadruple,interval_length,'quadruple')
    return single,double,triple,quadruple