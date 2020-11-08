# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:35:59 2020

@author: Лиля
"""

import os
import numpy as np
import pandas as pd

from scipy.stats import chisquare

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
                
    #результат нужно сохранить
    p_value.to_csv(os.path.join(results_saving_dir, 'p_value.csv'), sep='\t')

    logging.debug("P-value matrix\n%s", p_value)
    
    # logging.debug(msg=u'p-value matrix was created')
    
    return p_value

def primary_motifs(occurrences,background_n,p_value,acids,args,results_saving_dir):
    
    result=p_value[p_value<args.p_value/occurrences.size]*occurrences[occurrences>args.occurrences]
    result = pd.DataFrame(np.array([np.array(occurrences.index)] * occurrences.shape[1]).T, index=occurrences.index, columns=occurrences.columns)
    
    table = table_creator()
    
    for i in result[(p_value < args.p_value) & (occurrences > args.occurrences)].apply(lambda x: x.dropna().to_dict(), axis=1):
        if i:
            for pos, aa in i.items():
                if aa!='-':
                    n_motif = np.zeros(args.interval_length*2 + 1)
                    n_motif[pos+args.interval_length] = acids.index(aa) + 1
                    l_motif = ['.']*(args.interval_length*2 + 1)
                    l_motif[args.interval_length], l_motif[pos+args.interval_length]=args.modification_site.lower(), aa
    
                    motif = ''.join(l_motif).strip('.')
                    
                    table = table.append({'Letter motif': motif,'Number motif':n_motif,
                    'Observed':occurrences[pos][aa],'Expected':'-','Fasta':background_n[pos][aa],
                                        'p-value':p_value[pos][aa]},
                                            ignore_index = True)

    vector = np.array(table['Number motif'].tolist())
    logging.debug("Primary motifs:\n%s", table['Letter motif'])
    utils.saving_table(results_saving_dir,table,args.interval_length,'primary')
               
    return vector,table         



def counter(args, vector, acid_location, acid_number, dataset_info, previous_info, acids=utils.ACIDS_LIST):

    int_table, back_table = dataset_info 
    for i,ik in zip(acid_location, acid_number):
        int_table=int_table[(int_table[i-args.interval_length]==acids[ik-1])]
        back_table=back_table[(back_table[i-args.interval_length]==acids[ik-1])]         
        # print('acid',acids[ik-1])
    prev_observed, prev_fasta = previous_info    
    observed=len(int_table)
    fasta = len(back_table)
    expected=(fasta/prev_fasta)*prev_observed
    return observed,expected,fasta

    
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
           
def chi2_motifs(args, vector, acid_location, acid_number , dataset_info, previous_info):
    

    observed,expected, fasta = counter(args, vector, acid_location, acid_number, dataset_info, previous_info, acids=utils.ACIDS_LIST)
    if fasta==0:
        p_value = None
    else:
        chisq,p_value=chisquare([observed,previous_info[0]-observed],f_exp=[expected,previous_info[1]-expected])

    return observed, expected, fasta, p_value

def acids_n_l(motif):
    
    acid_location = [elem for elem in np.nonzero(motif)[0]]
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

def motif_counter(args,vector,motif,result,dataset_info, previous_info):
    acid_location, acid_number = acids_n_l(motif)
    observed, expected, fasta, p_value=chi2_motifs(args, vector, acid_location, acid_number , dataset_info, previous_info)
    if p_value!=None:
        motif_l=letter_motif(args,acid_location, acid_number,acids=utils.ACIDS_LIST)
        # print('motif_l',motif_l)
        result=result.append({'Letter motif':motif_l,'Number motif':motif,
                        'Observed':observed,'Expected':expected,'Fasta':fasta,
                                            'p-value':p_value},
                                                ignore_index = True)
    
    return result

def double_motifs(vector, dataset_info, results_saving_dir, args, single, acids=utils.ACIDS_LIST): 
    
    matrix = multiplying_matrix(args, vector, vector)                
    #создаем пустую табличку, в которую будем записывать результаты
    result = table_creator()
   
    for i in range(len(vector)):
        j=0
        while j<=i:
            elem=matrix[i,j]
            #нужны элементы матрицы с одними нулями
            if (elem.any())==False:
                motif=vector[i]+vector[j]
                prev_fg_occ = single['Observed'][i]
                prev_bg_occ = single[ 'Fasta'][i]
                previous_info = prev_fg_occ, prev_bg_occ
                result = motif_counter(args,vector,motif,result,dataset_info, previous_info)

            j+=1
            
    b=len(result['Observed'].values)
    table = (result.loc[result['p-value']<args.p_value/b][result['Observed']>=args.occurrences]).reset_index()        
    
    logging.debug("Double motifs:\n%s", table['Letter motif'])
    utils.saving_table(results_saving_dir,table,args.interval_length,'double')
                                                                     
    return table
   
def n_vectors(args, previous_motifs_table):
    
    vector=np.zeros((len(previous_motifs_table['Number motif'].values),args.interval_length*2+1))
    for i,ik in enumerate(previous_motifs_table['Number motif'].values):
        vector[i,:]=previous_motifs_table['Number motif'].values[i]
        
    return vector

def n_motifs_result(args, n_vector, vector, previous_table, matrix, result, dataset_info):
    for i in range(len(n_vector)):
        for j in range(len(vector)):
            elem=matrix[i,j]

            if (elem.any())==False:
                motif=vector[j]+n_vector[i]
            
                prev_fg_occ = previous_table['Observed'][i]
                prev_bg_occ = previous_table[ 'Fasta'][i]
                previous_info = prev_fg_occ, prev_bg_occ
              
                result = motif_counter(args,n_vector,motif,result,dataset_info,previous_info)
                
    b=len(result['Observed'].values)
    table = (result.loc[result['p-value']<args.p_value/b][result['Observed']>=args.occurrences]).reset_index()        
    return table           
    


def get_motifs(args, vector, previous_motifs_table, dataset_info, results_saving_dir, step):
    
    previous_vector = n_vectors(args, previous_motifs_table)
    matrix = multiplying_matrix(args, previous_vector, vector)
    result = table_creator()
    result = n_motifs_result(args, previous_vector, vector, previous_motifs_table, matrix, result, dataset_info)
    result.drop_duplicates(subset = ['Letter motif'], inplace = True)
    logging.info("%d - AA motifs:\n%s" % (step, result['Letter motif']))

    return result

def chi2_alg(fg_occ, bg_occ, dataset_info, args, results_saving_dir):
    
    p_values = p_value(fg_occ, bg_occ, args.interval_length, results_saving_dir)
    vector,single = primary_motifs(fg_occ, bg_occ ,p_values, utils.ACIDS_LIST, args, results_saving_dir)
    # dataset_info = fg_intervals, bg_intervals
    double=double_motifs(vector, dataset_info, results_saving_dir, args, single, acids=utils.ACIDS_LIST)   
    previous_motifs_table = double
    step = 3
    while True:
        new_motifs_table = get_motifs(args, vector, previous_motifs_table, dataset_info, results_saving_dir, step)

        if new_motifs_table.empty:
            break
        else:

            utils.saving_table(results_saving_dir, new_motifs_table, args.interval_length,' '.join([str(step),'AA motif']))                    
            previous_motifs_table = new_motifs_table
            step +=1
