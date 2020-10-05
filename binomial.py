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

def P_counter_bi(occurrences, background_n, args, results_saving_dir, acids=utils.ACIDS_LIST):

    
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
    
    
    logging.info(u'Binomial probability for each amino acid in matrix was counted')
    return P_binomial


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
    primary_motifs_probability=[]
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
                primary_motifs_probability.append(P_binomial[i][j])
    vector=np.array(primary_motifs_number)    
    table=pd.DataFrame({'Number motif':primary_motifs_number,'Letter motif':primary_motifs_letter,'Probability':primary_motifs_probability})
    print(table)
    print(vector)      
    utils.saving_table(results_saving_dir,table,args.interval_length,'primary')
    logging.info(str(len(table['Number motif'].values))+u' primary (one acid length) motifs were identificated')
   
    return vector,table 

def counter(args, table, acid_location, acid_number, acids = utils.ACIDS_LIST):
    for i, ik in zip(acid_location, acid_number):
        table = table[(table[i-args.interval_length]==acids[ik-1])]
    count = len(table)
    return count    

def l_motif(args,acid_location, acid_number,acids=utils.ACIDS_LIST):
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

def double_motifs_creator_bi(args, vector, single, intervals, background, acids = utils.ACIDS_LIST):
    
    table=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]), 'Observed': np.array([]),
                                                'P_AB':np.array([]),'P_A|B':np.array([])}) 
    
    b = np.tensordot(vector, vector.T, axes=0)
    l=len(vector)    
    matrix=np.zeros((l, l, args.interval_length * 2 + 1))
    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):      
            matrix[i,j,:] = np.diag(b[i,:,:,j])
    
    int_table = pd.DataFrame([list(i) for i in intervals], columns=range(-args.interval_length, args.interval_length + 1))
    back_table = pd.DataFrame([list(i) for i in background], columns=range(-args.interval_length, args.interval_length + 1))
    back_len = len(back_table)
    int_len = len(int_table)

    P_if = np.zeros((l, l, args.interval_length * 2 + 1))
    
    p_value = np.zeros((l, l, args.interval_length * 2 + 1))
    for i in range(l):
#        j=0
        for j in range(l):    

            elem=matrix[i,j]
        #нужны элементы матрицы с одними нулями
            if (elem.any())==False:

                motif = vector[i] + vector[j]
                acid_location = [elem for elem in np.nonzero(motif)[0]]
                acid_number = [int(motif[elem]) for elem in np.nonzero(motif)[0]]
                motif_l = l_motif(args, acid_location, acid_number,acids=utils.ACIDS_LIST)
            
                

                n_AB = counter(args, back_table, acid_location, acid_number, acids = utils.ACIDS_LIST)
                p_value = n_AB/back_len
                c = n_AB
                result = 0
                while c<=int_len:
                    result=result+binom.pmf(c,int_len,p_value,loc=0)
                    c+=1
                probability_j = single['Probability'][j]
                
                
                    
                P_if[i][j]=result/probability_j
                
                table=table.append({'Letter motif':motif_l,'Number motif':motif,'Observed':n_AB,
                                        'P_AB':result,'P_A|B':P_if[i][j]},
                                            ignore_index=True)
                

    print(table)
    return matrix
 
def triple_motifs_creator_bi(args, double_motifs):
    
#    table=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),'Observed': np.array([]),
#                                                'P_ABC':np.array([]),'P_A|BC':np.array([])})
    print(double_motifs)
    print(double_motifs['Observed'])
    b=len(double_motifs['Observed'].values)
    table=(double_motifs.loc[double_motifs['p-value']<args.p_value/b][double_motifs['Observed']>=args.occurrences]).reset_index()
    del table['index']
    
    vector=np.zeros((len(table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        vector[i,:]=table['Number motif'].values[i]
    print(vector)
    return vector                
#def triple_motifs_creator_bi(double_motifs, single_motifs, background, intervals, modification_site, interval_length, results_saving_dir, acids=utils.ACIDS_LIST):
#
#    indexes=[]
#    location=[]
#    probability=[]
#    occurrences=[]
#    backgrounds=[]
#    P_triple=[]
#    P_BC_triple=[]
#    #зафиксировали двубуквенную часть
#    for k in range(len(double_motifs['Location'].values)):
#        elem=double_motifs['Location'].values[k]
#        double_motif_1,double_motif_2=elem
#        double_motif_1_x,double_motif_1_y=double_motif_1
#        double_motif_2_x,double_motif_2_y=double_motif_2
#        
#        #теперь хотим добавлять однобуквенную 
#        #double_motif_AA=double_motifs['Location'].index[k]
#        for i in range(len(single_motifs['Location'].values)):
#            third_x,third_y=single_motifs['Location'].values[i]
#            #third_AA=acids[third_x]
#            #порверяем, что выбранная кислота не совпадает с предыдущими и нет накладок
#            if (third_y!=double_motif_1_y) and (third_y!=double_motif_2_y):
#                
#                    #считаем количество интервалов с тройным мотивом в background
#                    x=sum(1 for interval in background if ((interval[double_motif_1_y]==acids[double_motif_1_x]) 
#                                                           and (interval[double_motif_2_y]==acids[double_motif_2_x]) 
#                                                           and (interval[third_y]==acids[third_x])))
#
#                    #считаем p-value этого комплексного мотива
#                    p=x/len(background)
#
#                    #теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
#                    n=sum(1 for interval in intervals if ((interval[double_motif_1_y]==acids[double_motif_1_x]) 
#                                                           and (interval[double_motif_2_y]==acids[double_motif_2_x]) 
#                                                           and (interval[third_y]==acids[third_x])))
#
#                    #считаем Р вероятность по биномиальной формуле для тройного мотива
#                    P_ABC=0
#                    c=n
#                    while c<=len(intervals):
#                        P_ABC=P_ABC+binom.pmf(n,len(intervals),p,loc=0)
#                        c+=1
#                    #print('P_AB:',P_AB)
#
#                    #нашли вероятность P(ABC), теперь найдем условную вероятность
#                    P_BC=double_motifs['P_AB'].values[k]      
#                    #print('P_B:',P_B,motif_1_x,motif_1_y)
#                    P_A_BC=P_ABC/P_BC
#                    #print('P_A_B:',P_A_B)
#                
#                    #теперь разбираемся с названием мотива
#                    names=dict()
#                    names[double_motif_1_y]=acids[double_motif_1_x]
##                    print('names dict',double_motif_1_y,acids[double_motif_1_x])
#                    names[double_motif_2_y]=acids[double_motif_2_x]
##                    print('names dict',double_motif_2_y,acids[double_motif_2_x])
#                    names[third_y]=acids[third_x]
##                    print('names',third_y,acids[third_x])
#                    names[interval_length]=modification_site.lower()
#                    list_names=list(names.keys())
##                    print('list names 0',list_names)
#                    list_names.sort()
##                    print('list names',list_names)
#                    
#                    name=[names[list_names[0]],'.'*(list_names[1]-list_names[0]-1),names[list_names[1]],
#                          '.'*(list_names[2]-list_names[1]-1),names[list_names[2]],'.'*(list_names[3]-list_names[2]-1),
#                          names[list_names[3]]]
#                    motif=''.join(name)
#                          
#                    loc=((acids.index(names[list_names[0]].upper()),list_names[0]), (acids.index(names[list_names[1]].upper()),list_names[1]),
#                          (acids.index(names[list_names[2]].upper()),list_names[2]),(acids.index(names[list_names[3]].upper()),list_names[3]))
#                          
#                    probability.append(P_A_BC)
#
#                    occurrences.append(n)
#                    
#                    backgrounds.append(x)
#
#                    P_triple.append(P_ABC)
#                    
#                    P_BC_triple.append(P_BC)
#                          
#                    location.append(loc)
#                          
#                    indexes.append(motif)      
#                    
#                    
#    motifs=pd.DataFrame({'Location': location, 'Probability': probability , 'Occurrences' : occurrences, 'Backgrounds':backgrounds, 'P_ABC' : P_triple, 'P_BC' : P_BC_triple}, index=indexes)
#    sorted_motifs=motifs.copy()
#    sorted_motifs.sort_values('Probability',ascending=True,inplace=True)
#    triple_motifs=sorted_motifs
#    
#    count=[]
#    for i in range(len(triple_motifs['Probability'].values)):
#        if ((triple_motifs['Probability'].values[i])<0.05/len(triple_motifs['Probability'].values) and (triple_motifs['Occurrences'].values[i])>10):
#            count.append(1)
#        else:
#            count.append(0)
#    triple_motifs['Count']=count
#    
#    triple_motifs_selection=triple_motifs.where(triple_motifs['Count']==1).dropna()
#    triple_motifs_selection_copy=triple_motifs_selection.copy()
#    del triple_motifs_selection_copy['Count']
#       
#    if probability==[]:
#        triple_motifs_selection_copy=None
#    
##     if triple_motifs_selection_copy is not None:
## #         volcano_plot_for_motifs(triple_motifs_selection_copy,path_results,name='triple_motifs')    
#    logging.info(str(len(triple_motifs_selection_copy['Occurrences'].values))+u' triple (three acid length) motifs were identificated')    
#    return triple_motifs_selection_copy

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
    double = double_motifs_creator_bi(args, vector, single, intervals, background, acids = utils.ACIDS_LIST)
    triple = triple_motifs_creator_bi(args, double)
#    utils.saving_table(results_saving_dir,single,interval_length,'single')
#    double=double_motifs_creator_bi(args, vector, single, intervals, background, P_binomial, results_saving_dir, acids=utils.ACIDS_LIST)
#    if double is not None:
#        utils.saving_table(results_saving_dir,double,interval_length,'double')
#        triple=triple_motifs_creator_bi(double,single,background,intervals,modification_site,interval_length,results_saving_dir, acids=utils.ACIDS_LIST)
#    if triple is not None:
#        utils.saving_table(results_saving_dir,triple,interval_length,'triple')
#        quadruple=quadruple_motifs_creator_bi(triple,single,background,intervals,modification_site,interval_length,results_saving_dir, acids=utils.ACIDS_LIST)
#    if quadruple is not None:
#        utils.saving_table(results_saving_dir,quadruple,interval_length,'quadruple')
#    return single,double,triple,quadruple
    return vector,single, double, triple