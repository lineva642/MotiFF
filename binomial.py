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


def P_counter_bi(occurrences, background_n, args, results_saving_dir):
    fg_size = occurrences[0].sum()
    bg_prob = background_n / background_n[0].sum()#х
    binom_prob = pd.DataFrame(binom.sf(occurrences - 1, fg_size, bg_prob), columns=bg_prob.columns, index=bg_prob.index)
    utils.saving_table(results_saving_dir, binom_prob, args.interval_length, 'P_binomial_matrix')            
    logging.info(u'Binomial probability for each amino acid in matrix was counted')
    logging.debug("P_binomial matrix:\n%s", binom_prob)
    return binom_prob


# def letter_motif(args,indexes,acids=utils.ACIDS_LIST):
#     position=dict()
#     for elem in indexes:
#         position[elem]=acids[indexes[elem]-1]
#     position[args.interval_length]=(args.modification_site).lower()    

#     keys=list(position.keys())
#     keys.sort()
#     motif=position[keys[0]]+'.'*(keys[1]-keys[0]-1)
#     i=1
#     while i<len(keys)-1:
#         motif=''.join([motif,position[keys[i]],'.'*(keys[i+1]-keys[i]-1)])
#         i+=1
#     motif=''.join([motif,position[keys[len(keys)-1]]])
#     return motif

def get_letter_motif(args, acid_location, acid_number, acids=utils.ACIDS_LIST):
    motif = []
    position=dict()
    for i,ik in zip(acid_location, acid_number):
        position[i]=acids[ik-1]
    position[args.interval_length]= args.modification_site.lower()    

    keys=list(position.keys())
    keys.sort()
    motif=position[keys[0]]+'.'*(keys[1]-keys[0]-1)
    i=1
    while i<len(keys)-1:
        motif=''.join([motif,position[keys[i]],'.'*(keys[i+1]-keys[i]-1)])
        i+=1
    motif=''.join([motif,position[keys[len(keys)-1]]])
    return motif       

def single_motifs_creator_bi(binom_prob, occurrences, intervals, args, results_saving_dir, acids=utils.ACIDS_LIST):
    binom_prob = binom_prob[binom_prob < args.p_value / occurrences.size] # check this
    logging.debug('Single motif binomial probabilities:\n%s', binom_prob)
    occurrences = occurrences[occurrences > args.occurrences]
    result = binom_prob* occurrences
    logging.debug('Single motif result:\n%s', result)
    primary_motifs_number = []
    primary_motifs_letter = []
    primary_motifs_probability = []
    for i in result.columns:
        for j in result.index:
            if (math.isnan(result[i][j])):
                continue
            else:
                #n_motif-number motif, l_motif-letter motif
                n_motif=np.zeros(args.interval_length*2+1)
                n_motif[i+args.interval_length]=acids.index(j)+1
                
                acid_location, acid_number = [i+args.interval_length], [acids.index(j) + 1] 

                # indexes={(i+args.interval_length):(acids.index(j)+1)}
                l_motif= get_letter_motif(args,acid_location, acid_number,acids=utils.ACIDS_LIST)
                    
                primary_motifs_letter.append(l_motif)
                primary_motifs_number.append(n_motif)
                primary_motifs_probability.append(binom_prob[i][j])
    vector=np.array(primary_motifs_number)    
    table=pd.DataFrame({'Number motif':primary_motifs_number,'Letter motif':primary_motifs_letter,'Probability':primary_motifs_probability})
    # print(table)
    # print(vector)  
    logging.debug('Table %s', table)
    utils.saving_table(results_saving_dir,table,args.interval_length,'primary')
    logging.info(str(len(table['Number motif'].values))+u' primary (one acid length) motifs were identificated')
    print('vector', vector)
    return vector,table 

def counter(args, table, acid_location, acid_number, acids = utils.ACIDS_LIST):
    for i, ik in zip(acid_location, acid_number):
        table = table[(table[i-args.interval_length]==acids[ik-1])]
    count = len(table)
    return count    


def double_motifs_creator_bi(args, vector, single, intervals, background, acids = utils.ACIDS_LIST):
    
    table=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]), 'Observed': np.array([]),
                                                'P_AB':np.array([]),'P_A|B':np.array([])}) 
    #составляем всевозможные пары двойных мотивов из полученных первичных, нулевые элементы-нужные нам мотивы
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

    
    
    for i in range(l):
#        j=0
        for j in range(l):    

            elem=matrix[i,j]
        #нужны элементы матрицы с одними нулями
            if (elem.any())==False:
                # print('double_motif',motif)
                motif = vector[i] + vector[j]
                # print('double_motif',motif)
                #хотим восстановить буквенный вид мотива
                acid_location = [elem for elem in np.nonzero(motif)[0]]
                acid_number = [int(motif[elem]) for elem in np.nonzero(motif)[0]]
                motif_l = letter_motif(args, acid_location, acid_number,acids=utils.ACIDS_LIST)
            
                

                n_AB = counter(args, back_table, acid_location, acid_number, acids = utils.ACIDS_LIST)
                p_value = n_AB/back_len
                c = n_AB
                result = 0
                while c<=int_len:
                    result=result+binom.pmf(c,int_len,p_value,loc=0)
                    c+=1
                probability_j = single['Probability'][j]
                
                
                    
                P_if=result/probability_j
                
                table=table.append({'Letter motif':motif_l,'Number motif':motif,'Observed':n_AB,
                                        'P_AB':result,'P_A|B':P_if},
                                            ignore_index=True)
                

    # print(table)
    return table
 
def triple_motifs_creator_bi(args, vector, double_motifs, intervals, background):
    
#    table=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]),'Observed': np.array([]),
#                                                'P_ABC':np.array([]),'P_A|BC':np.array([])})
    result_table=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]), 'Observed': np.array([]),
                                                'P_ABC':np.array([]),'P_A|BC':np.array([])}) 
    
    b=len(double_motifs['Observed'].values)
    table=(double_motifs.loc[double_motifs['P_A|B']<args.p_value/b][double_motifs['Observed']>=args.occurrences]).reset_index()
    del table['index']
    
    # print('!!!',table)

    
    double_vector=np.zeros((len(table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        double_vector[i,:]=table['Number motif'].values[i]
    # print('!!!',double_vector)
    
    b=np.tensordot(double_vector,vector.T,axes=0)
    matrix=np.zeros((len(double_vector),len(vector),args.interval_length*2+1))

    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j])
    # print('!!!',matrix)

    int_table = pd.DataFrame([list(i) for i in intervals], columns=range(-args.interval_length, args.interval_length + 1))
    back_table = pd.DataFrame([list(i) for i in background], columns=range(-args.interval_length, args.interval_length + 1))
    back_len = len(back_table)
    int_len = len(int_table)

    
    
    for i in range(len(double_vector)):
#        j=0
        for j in range(len(vector)):    

            elem=matrix[i,j]
        #нужны элементы матрицы с одними нулями
            if (elem.any())==False:

                motif = double_vector[i] + vector[j]
                #хотим восстановить буквенный вид мотива
                acid_location = [elem for elem in np.nonzero(motif)[0]]
                acid_number = [int(motif[elem]) for elem in np.nonzero(motif)[0]]
                motif_l = letter_motif(args, acid_location, acid_number,acids=utils.ACIDS_LIST)
            
                

                n_ABC = counter(args, back_table, acid_location, acid_number, acids = utils.ACIDS_LIST)
                p_value = n_ABC/back_len
                c = n_ABC
                result = 0
                while c<=int_len:
                    result=result+binom.pmf(c,int_len,p_value,loc=0)
                    c+=1
                probability_i = table['P_AB'][i]
                
                
                    
                P_if=result/probability_i
                
                result_table=result_table.append({'Letter motif':motif_l,'Number motif':motif,'Observed':n_ABC,
                                        'P_ABC':result,'P_A|BC':P_if},
                                            ignore_index=True)
    # print('!!! result',result_table)            
        
    return result_table   

def quadruple_motifs_creator_bi(args, vector, triple_motifs, intervals, background):

    result_table=pd.DataFrame({'Letter motif':np.array([]),'Number motif':np.array([]), 'Observed': np.array([]),
                                                'P_ABCD':np.array([]),'P_A|BCD':np.array([])}) 
    
    b=len(triple_motifs['Observed'].values)
    table=(triple_motifs.loc[triple_motifs['P_A|BC']<args.p_value/b][triple_motifs['Observed']>=args.occurrences]).reset_index()
    del table['index']
    
    # print('!!!!',table)
    
    
    triple_vector=np.zeros((len(table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(table['Number motif'].values):
        triple_vector[i,:]=table['Number motif'].values[i]
    # print('!!!!',triple_vector)
    
    b=np.tensordot(triple_vector,vector.T,axes=0)
    matrix=np.zeros((len(triple_vector),len(vector),args.interval_length*2+1))

    for i,ik in enumerate(b):
        for j,jk in enumerate(b[0][0][0]):
            matrix[i,j,:]=np.diag(b[i,:,:,j])
    # print('!!!!',matrix)

    int_table = pd.DataFrame([list(i) for i in intervals], columns=range(-args.interval_length, args.interval_length + 1))
    back_table = pd.DataFrame([list(i) for i in background], columns=range(-args.interval_length, args.interval_length + 1))
    back_len = len(back_table)
    int_len = len(int_table)

    
    
    for i in range(len(triple_vector)):
#        j=0
        for j in range(len(vector)):    

            elem=matrix[i,j]
        #нужны элементы матрицы с одними нулями
            if (elem.any())==False:
                # print('quadruple_motif',motif)
                motif = triple_vector[i] + vector[j]
                print('quadruple_motif',motif)
                #хотим восстановить буквенный вид мотива
                acid_location = [elem for elem in np.nonzero(motif)[0]]
                acid_number = [int(motif[elem]) for elem in np.nonzero(motif)[0]]
                motif_l = letter_motif(args, acid_location, acid_number,acids=utils.ACIDS_LIST)
            
                

                n_ABCD = counter(args, back_table, acid_location, acid_number, acids = utils.ACIDS_LIST)
                p_value = n_ABCD/back_len
                c = n_ABCD
                result = 0
                while c<=int_len:
                    result=result+binom.pmf(c,int_len,p_value,loc=0)
                    c+=1
                probability_i = table['P_ABC'][i]
                
                
                    
                P_if=result/probability_i
                
                result_table=result_table.append({'Letter motif':motif_l,'Number motif':motif,'Observed':n_ABCD,
                                        'P_ABCD':result,'P_A|BCD':P_if},
                                            ignore_index=True)
    # print('!!!! result',result_table)            
        
    return triple_vector    
    

def motifs_bi(binom_prob, occurrences, idPeptides, background, args, results_saving_dir):

    intervals = idPeptides['fasta_match']
    vector, single = single_motifs_creator_bi(binom_prob, occurrences, intervals, args,  results_saving_dir, acids=utils.ACIDS_LIST)
    utils.saving_table(results_saving_dir,single,args.interval_length,'single')
    double = double_motifs_creator_bi(args, vector, single, intervals, background, acids = utils.ACIDS_LIST)
    if double is not None:
        b=len(double['Observed'].values)
        result_double = (double.loc[double['P_A|B']<args.p_value/b][double['Observed']>=args.occurrences]).reset_index()
        utils.saving_table(results_saving_dir,result_double,args.interval_length,'double')
        triple = triple_motifs_creator_bi(args, vector, double, intervals, background)
        if triple is not None:
            b=len(triple['Observed'].values)
            result_triple = (triple.loc[triple['P_A|BC']<args.p_value/b][triple['Observed']>=args.occurrences]).reset_index()
            utils.saving_table(results_saving_dir,result_triple,args.interval_length,'triple')
            quadruple = quadruple_motifs_creator_bi(args, vector, triple, intervals, background)
            if quadruple is not None:
                b=len(quadruple['Observed'].values)
                result_quadruple = (quadruple.loc[quadruple['P_A|BC']<args.p_value/b][quadruple['Observed']>=args.occurrences]).reset_index()
                utils.saving_table(results_saving_dir,result_quadruple,args.interval_length,'quadruple')
            else:
                quadruple = None
        else:
            triple = None
    else:
        double = None           
    return single, double, triple, quadruple