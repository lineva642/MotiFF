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

#def p_value_bi(modification_site,interval_length,background,background_n, acids=utils.ACIDS_LIST):
#    background_n_1=[]
#    for i in range(len(acids)):
#        if acids[i]==modification_site:
#            number=i
#    for i in range(len(acids)):
#        a=[0]*(interval_length*2+1)
#        background_n_1.append(a)
#    #составили матрицу вероятностей
#    for i in range(len(acids)):
#        for k in range(interval_length*2+1):
#            if (k==interval_length) and (i==number):
#                background_n_1[i][k]=1
#            else:    
#                background_n_1[i][k]=(background_n[i][k])/len(background)
#    logging.debug(msg=u'p for each background occurrence matrix element was counted')            
#    return background_n_1

#def P_matrix_bi(interval_length,n,N,background_n_1,results_saving_dir, acids=utils.ACIDS_LIST):
#    P=[]
#    path=os.path.join(results_saving_dir,'P.txt')
#    saving=open(path,'w')
#    for i in range(len(acids)):
#        a=[0]*(interval_length*2+1)
#        P.append(a)
#    for i in range(len(acids)):
#        for k in range(interval_length*2+1):
##            print(i,k)
#            #формула-binom.pmf(k, n, p, loc=0)
#            result=0
#            c=n[i][k]
#            while c<=N:
#                result=result+binom.pmf(c,N,background_n_1[i][k],loc=0)
#                c+=1
#            P[i][k]=result
#            saving.write(str(result)+' ')
#        saving.write('\n')
#    saving.close()
#    logging.debug(msg=u'Binomal probability for each occurrence matrix element was counted')    
#    return P

def P_counter_bi(occurrences, background_n, args, results_saving_dir, acids=utils.ACIDS_LIST):
#    print('Probability is being counted')
#    n = utils.n_calculator(intervals, interval_length, results_saving_dir,acids=utils.ACIDS_LIST)
#    N_calculator = lambda intervals:len(intervals)
#    N = N_calculator(intervals)
#    background_n = utils.background_n_matrix(interval_length,background,results_saving_dir, acids=utils.ACIDS_LIST)
    
    p_value = background_n/(background_n[background_n.columns[0]].sum())
    P_binomial=pd.DataFrame().reindex_like(p_value)
    for i in range( -args.interval_length, args.interval_length + 1):
        for acid in acids:
            result=0
            c=occurrences[i][acid]
            while c<=(occurrences[occurrences.columns[0]].sum()):
                result=result+binom.pmf(c,occurrences[occurrences.columns[0]].sum(),p_value[i][acid],loc=0)
                c+=1
            P_binomial[i][acid]=result
    utils.saving_table(results_saving_dir,P_binomial,args.interval_length,'P_binomial_matrix')            
    print(P_binomial)        
    
    
#    background_n_1=p_value_bi(modification_site,interval_length,background,background_n, acids=utils.ACIDS_LIST)
#    P = P_matrix_bi(interval_length,n,N,background_n_1,results_saving_dir, acids=utils.ACIDS_LIST)
#    print('Probability is being counted')
    logging.info(u'Binomial probability for each amino acid in matrix was counted')
    return P_binomial

#рассчитаем occurrences
#def occurrences_counter_bi(intervals, interval_length, modification_site, results_saving_dir, acids=utils.ACIDS_LIST):
#    n=utils.n_calculator(intervals, interval_length, results_saving_dir, acids=utils.ACIDS_LIST)
##    N=N_calculator(intervals)
#    path=os.path.join(results_saving_dir,'occurrences.txt')
#    saving=open(path,'w')
#    occurrences=[]
#    for i in range(len(acids)):
#        a=[0]*(interval_length*2+1)
#        occurrences.append(a)
#    for i in range(len(acids)):
#        for k in range(interval_length*2+1):
#            if k!=interval_length:
#                occurrences[i][k]=n[i][k]
#            saving.write(str(occurrences[i][k])+' ')
#        saving.write('\n')
#    saving.close()             
#    logging.debug(u'Occurrences for each amino acid in matrix was counted')            
#    return occurrences 

#отбор по вероятности
#def P_validation_bi(interval_length, P, acids=utils.ACIDS_LIST):
#    P1=[]
#    for i in range(len(acids)):
#        a=[0]*(interval_length*2+1)
#        P1.append(a)
#    for i in range(len(acids)):
#        for k in range(interval_length*2+1):
#            if P[i][k]<0.05/(21*13):
#                continue
#            else:
#                P1[i][k]+=1         
#    return P1
#
##отбор по встречаемости более 10
#def occurrences_validation_bi(interval_length, occurrences, acids=utils.ACIDS_LIST):
#    occurrences_1=[]
#    for i in range(len(acids)):
#        a=[0]*(interval_length*2+1)
#        occurrences_1.append(a)
#    for i in range(len(acids)):
#        for k in range(interval_length*2+1):
#            if occurrences[i][k]>=10:
#                occurrences_1[i][k]=occurrences[i][k]
#    return occurrences_1
#
##соединим валидацию по occurances и p-values
##нужно совместить occurances_1 & P1
#def final_validation_bi(interval_length, occurrences, P, acids=utils.ACIDS_LIST):
#    occurrences_1=occurrences_validation_bi(interval_length, occurrences, acids=utils.ACIDS_LIST)
#    P1=P_validation_bi(interval_length, P, acids=utils.ACIDS_LIST)
#    for i in range(len(acids)):
#        for k in range(interval_length*2+1):
#            if occurrences_1[i][k]>=10:
#                continue
#            else:
#                P1[i][k]=1
#    (u'Binomial probabilities were counted and selection with p=0.05/Bonferroni correction and occurrences>10 was performed')            
#    return P1 

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

def single_motifs_creator_bi(args, P_binomial, occurrences, intervals, results_saving_dir, acids=utils.ACIDS_LIST):
    P_binomial=P_binomial[P_binomial<args.p_value/occurrences.size]
    occurrences=occurrences[occurrences>args.occurrences]
    result=P_binomial*occurrences
    
    primary_motifs_number=[]
    primary_motifs_letter=[]
    for i in result.columns:
        for j in result.index:
            if (math.isnan(result[i][j])):
                continue
            else:
                #n_motif-number motif, l_motif-letter motif
                n_motif=np.zeros(args.interval_length*2+1)
                n_motif[i+args.interval_length]=acids.index(j)+1

                indexes={(i+args.interval_length):(acids.index(j)+1)}
                l_motif=letter_motif(args,indexes,acids=utils.ACIDS_LIST)
                    
                primary_motifs_letter.append(l_motif)
                primary_motifs_number.append(n_motif)
    vector=np.array(primary_motifs_number)    
    table=pd.DataFrame({'Number motif':primary_motifs_number,'Letter motif':primary_motifs_letter})
    print(table)
    print(vector)      
    utils.saving_table(results_saving_dir,table,args.interval_length,'primary')
    logging.info(str(len(table['Number motif'].values))+u' primary (one acid length) motifs were identificated')
   
    return vector,table 


def double_motifs_creator_bi(args, vector, intervals, background, P_binomial, results_saving_dir, acids=utils.ACIDS_LIST):
    
    b = np.tensordot(vector, vector.T, axes=0)
    l=len(vector)    
    matrix=np.zeros((l, l, args.interval_length * 2 + 1))
    for i,ik in enumerate(b):    
        j= 0    
        while j<= i:        
            matrix[i,j,:] = np.diag(b[i,:,:,j])
            j+=1 

    print('matrix',matrix)               
    #создаем пустую табличку, в которую будем записывать результаты
    result=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),
                                                'Observed':np.array([]),'Expected':np.array([]),
                                                                    'p-value':np.array([])})
#    intervals=(idPeptides['fasta_match']).sum()
    
    int_len= len(intervals)
    back_len= len(background)
    
    for i in range(l):
        j=0
        while j<=i:
            elem=matrix[i,j]
            #нужны элементы матрицы с одними нулями
            if (elem.any())==False:
                motif=vector[i]+vector[j]
                print(motif)
            j+=1    
    return result
    
    
#    #для каждого претендента нужно построить двойной комплексный мотив и оценить вероятность такого мотива
#    indexes=[]
#    location=[]
#    probability=[]
#    occurrences=[]
#    backgrounds=[]
#    P_double=[]
#    P_B_double=[]
#
#        
#    #закрепили одну аминокислоту
#    for i in range(len(single_motifs_creator)):
#        #выбрали вторую аминокислоту
#        #print(i)
#        for j in range(len(single_motifs_creator)):
#            if j!=i:
#                #print(j)
#                #выбрали наш сложный мотив как две АА
#                motif_1_x,motif_1_y=single_motifs_creator['Location'].values[i]
#                motif_2_x,motif_2_y=single_motifs_creator['Location'].values[j]
#
#                #рассматриваем кислоты, которые стоят в интервале на разных позициях
#                if motif_1_y!=motif_2_y:
#                    
#                    #считаем количество интервалов с данным комплексным мотивом в background
#                    x=sum(1 for interval in background if ((interval[motif_1_y]==acids[motif_1_x]) and (interval[motif_2_y]==acids[motif_2_x])))
#
#                    #считаем p-value этого комплексного мотива
#                    p=x/len(background)
#
#                    #теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
#                    n=sum(1 for interval in intervals if ((interval[motif_1_y]==acids[motif_1_x]) and (interval[motif_2_y]==acids[motif_2_x])))
#
#                    #считаем Р вероятность по биномиальной формуле для комплексного мотива
#                    P_AB=0
#                    c=n
#                    while c<=len(intervals):
#                        P_AB=P_AB+binom.pmf(n,len(intervals),p,loc=0)
#                        c+=1
#                    #print('P_AB:',P_AB)
#
#                    #нашли вероятность P(AB), теперь найдем условную вероятность
#                    P_B=P[motif_1_x][motif_1_y]
#                    #print('P_B:',P_B,motif_1_x,motif_1_y)
#                    P_A_B=P_AB/P_B
#                    #print('P_A_B:',P_A_B)
#
#                    
#                    #создадим dataframe
#                    if (motif_1_y<motif_2_y) and (motif_1_y<interval_length) and (motif_2_y>interval_length):
#                        s = [acids[motif_1_x], '.'*(interval_length-motif_1_y-1), modification_site.lower(),'.'*(motif_2_y-interval_length-1),acids[motif_2_x]]
#                        index=''.join(s)
#                        #s=motif_1_acid+'.'*(interval_length-motif_1_y-1)+modification_site.lower()+'.'*(motif_2_y-interval_length-1)+motif_2_acid
#                        indexes.append(index)
#                    elif (motif_1_y<motif_2_y) and (motif_1_y<interval_length) and (motif_2_y<interval_length):
#                        s = [acids[motif_1_x], '.'*(motif_2_y-motif_1_y-1),acids[motif_2_x],'.'*(interval_length-motif_2_y-1),modification_site.lower()]
#                        index=''.join(s)
#                        #s=motif_1_acid+'.'*(motif_2_y-motif_1_y-1)+motif_2_acid+'.'*(interval_length-motif_2_y-1)+modification_site.lower()
#                        indexes.append(index)
#                    elif (motif_1_y<motif_2_y) and (motif_1_y>interval_length):
#                        s = [modification_site.lower(), '.'*(motif_1_y-interval_length-1),acids[motif_1_x],'.'*(motif_2_y-motif_1_y-1),acids[motif_2_x]]
#                        index=''.join(s)
#                        #s=modification_site.lower()+'.'*(motif_1_y-interval_length-1)+motif_1_acid+'.'*(motif_2_y-motif_1_y-1)+motif_2_acid
#                        indexes.append(index)
#                    elif (motif_2_y<motif_1_y) and (motif_2_y<interval_length) and (motif_1_y>interval_length):
#                        s = [acids[motif_2_x], '.'*(interval_length-motif_2_y-1),modification_site.lower(),'.'*(motif_1_y-interval_length-1),acids[motif_1_x]]
#                        index=''.join(s)
#                        #s=motif_2_acid+'.'*(interval_length-motif_2_y-1)+modification_site.lower()+'.'*(motif_1_y-interval_length-1)+motif_1_acid
#                        indexes.append(index)
#                    elif (motif_2_y<motif_1_y) and (motif_2_y<interval_length) and (motif_1_y<interval_length):
#                        s = [acids[motif_2_x], '.'*(motif_1_y-motif_2_y-1),acids[motif_1_x],'.'*(interval_length-motif_1_y-1),modification_site.lower()]
#                        index=''.join(s)                        
#                        #s=motif_2_acid+'.'*(motif_1_y-motif_2_y-1)+motif_1_acid+'.'*(interval_length-motif_1_y-1)+modification_site.lower()
#                        indexes.append(index)
#                    elif (motif_2_y<motif_1_y) and (motif_2_y>interval_length):
#                        s = [modification_site.lower(), '.'*(motif_2_y-interval_length-1),acids[motif_2_x],'.'*(motif_1_y-motif_2_y-1),acids[motif_1_x]]
#                        index=''.join(s)                            
#                        #s=modification_site.lower()+'.'*(motif_2_y-interval_length-1)+motif_2_acid+'.'*(motif_1_y-motif_2_y-1)+motif_1_acid
#                        indexes.append(index)
#
#                    if motif_1_y<motif_2_y:
#                        location.append(((motif_1_x,motif_1_y),(motif_2_x,motif_2_y)))
#                    else:
#                        location.append(((motif_2_x,motif_2_y),(motif_1_x,motif_1_y)))
#                    #location.append((motif_2_x,motif_2_y))
#
#                    probability.append(P_A_B)
#
#                    occurrences.append(n)
#                    
#                    backgrounds.append(x)
#
#                    P_double.append(P_AB)
#                    
#                    P_B_double.append(P_B)
#                    
#                    
#    motifs=pd.DataFrame({'Location': location, 'Probability': probability , 'Occurrences' : occurrences, 'Backgrounds': backgrounds, 'P_AB' : P_double, 'P_B' : P_B_double}, index=indexes)
#    sorted_motifs=motifs.copy()
#    sorted_motifs.sort_values('Probability',ascending=True,inplace=True)
#    double_motifs=sorted_motifs
#    
#    count=[]
#    for i in range(len(double_motifs['Probability'].values)):
#        if ((double_motifs['Probability'].values[i])<0.05/len(double_motifs['Probability'].values) and (double_motifs['Occurrences'].values[i])>10):
#            count.append(1)
#        else:
#            count.append(0)
#    double_motifs['Count']=count
#    
#    double_motifs_selection=double_motifs.where(double_motifs['Count']==1).dropna()
#    double_motifs_selection_copy=double_motifs_selection.copy()
#    del double_motifs_selection_copy['Count']
#    
#    
##    print('double_motifs',double_motifs_selection_copy)
#    
#    if probability==[]:
#        double_motifs_selection_copy=None
#    
##     if double_motifs_selection_copy is not None:
##         volcano_plot_for_motifs(double_motifs_selection_copy,path_results,name='double_motifs')    
#        
#    logging.info(str(len(double_motifs_selection_copy['Occurrences'].values))+u' double (two acid length) motifs were identificated')
#    return double_motifs_selection_copy

def triple_motifs_creator_bi(double_motifs, single_motifs, background, intervals, modification_site, interval_length, results_saving_dir, acids=utils.ACIDS_LIST):

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

def quadruple_motifs_creator_bi(triple_motifs, single_motifs, background, intervals, modification_site, interval_length, results_saving_dir, acids=utils.ACIDS_LIST):

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

def motifs_bi(args, P_binomial, occurrences, idPeptides, background, results_saving_dir, acids=utils.ACIDS_LIST):
    #primary=primary_motifs_creator(acids,P_final,interval_length,modification_site)
    intervals=(idPeptides['fasta_match']).sum()
    vector,single=single_motifs_creator_bi(args, P_binomial, occurrences, intervals, results_saving_dir, acids=utils.ACIDS_LIST)
#    utils.saving_table(results_saving_dir,single,interval_length,'single')
    double=double_motifs_creator_bi(args, vector, idPeptides, background, P_binomial, results_saving_dir, acids=utils.ACIDS_LIST)
#    if double is not None:
#        utils.saving_table(results_saving_dir,double,interval_length,'double')
#        triple=triple_motifs_creator_bi(double,single,background,intervals,modification_site,interval_length,results_saving_dir, acids=utils.ACIDS_LIST)
#    if triple is not None:
#        utils.saving_table(results_saving_dir,triple,interval_length,'triple')
#        quadruple=quadruple_motifs_creator_bi(triple,single,background,intervals,modification_site,interval_length,results_saving_dir, acids=utils.ACIDS_LIST)
#    if quadruple is not None:
#        utils.saving_table(results_saving_dir,quadruple,interval_length,'quadruple')
#    return single,double,triple,quadruple
    return vector,single,double