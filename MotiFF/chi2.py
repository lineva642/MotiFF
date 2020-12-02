# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 23:35:59 2020

@author: Лиля
"""

import os
import numpy as np
import pandas as pd
import logging
import utils

from scipy.stats import chisquare

##АЛГОРИТМ_chi2

      
# считаем p-value
def calculate_p_value(fg_occ, bg_occ, fg_size, bg_size, interval_length):
    norm_expected = bg_occ / bg_size * fg_size
    print(f'{fg_size=}, {bg_size=}')
    p_value = fg_occ.combine(norm_expected, lambda fg, bg: chisquare([fg, fg_size - fg], [bg, fg_size - bg])[1])
    logging.debug("P-value matrix\n%s", p_value)
    return p_value


def primary_motifs(fg_occ, bg_occ, p_value, args, acids=utils.ACIDS_LIST):
    
    # HERE is multiple comparison should be counted
    # result = p_value[p_value < args.p_value] * fg_occ[fg_occ > args.occurrences] #p_value[p_value < args.p_value / occurrences.size] * occurrences[occurrences > args.occurrences] 
    # print(f"{result=}")
    result = pd.DataFrame(np.array([np.array(fg_occ.index)] * fg_occ.shape[1]).T, index=fg_occ.index, columns=fg_occ.columns)
    print(f"{result=}")
    table = table_creator()
    
    for i in result[(p_value < args.p_value) & (fg_occ > args.occurrences)].apply(lambda x: x.dropna().to_dict(), axis=1):
        if i:
            for pos, aa in i.items():
                if aa != '-':
                    n_motif = np.zeros(args.interval_length * 2 + 1)
                    n_motif[pos + args.interval_length] = acids.index(aa) + 1
                    l_motif = ['.'] * (args.interval_length * 2 + 1)
                    l_motif[args.interval_length] = args.modification_site.lower()
                    l_motif[pos + args.interval_length] = aa
                    motif = ''.join(l_motif).strip('.')
                    table = table.append({'motif': motif,
                                          'Number motif': n_motif,
                                          'p_value': p_value[pos][aa],
                                          'fg_matches': fg_occ[pos][aa],
                                          'fg_size': fg_occ[args.interval_length].sum(), 
                                          'bg_size': bg_occ[args.interval_length].sum(),
                                          'bg_matches': bg_occ[pos][aa]},
                                         ignore_index=True)

    vector = np.array(table['Number motif'].tolist())
    logging.debug("Primary motifs:\n%s", table['motif'])  
    print(f"{vector=} \n {table=}")       
    return vector, table         


def double_motifs(vector, dataset_info, args, single, acids=utils.ACIDS_LIST): 
    
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
                prev_fg_occ = single['fg_matches'][i]
                prev_bg_occ = single['bg_matches'][i]
                previous_info = prev_fg_occ, prev_bg_occ
                result = motif_counter(args, vector, motif, result, dataset_info, previous_info)
            j+=1
    b=len(result['fg_matches'].values)
    table = (result.loc[ result['p_value'] < args.p_value / b][result['fg_matches'] >= args.occurrences]).reset_index()        
    logging.debug("Double motifs:\n%s", table['motif'])
    print("Table 2 \n{table}")
    return table


def counter(args, vector, acid_location, acid_number, dataset_info, previous_info, acids=utils.ACIDS_LIST):
    int_table, back_table = dataset_info 
    for i,ik in zip(acid_location, acid_number):
        int_table=int_table[(int_table[i-args.interval_length]==acids[ik-1])]
        back_table=back_table[(back_table[i-args.interval_length]==acids[ik-1])]         
    prev_observed, prev_fasta = previous_info    
    observed=len(int_table)
    fasta = len(back_table)
    expected=(fasta/prev_fasta)*prev_observed
    return observed,expected,fasta

    
def letter_motif(args, acid_location, acid_number, acids=utils.ACIDS_LIST):
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
    observed, expected, fasta = counter(args, vector, acid_location, acid_number, dataset_info, previous_info, acids=utils.ACIDS_LIST)
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
    result = pd.DataFrame({'motif':np.array([]),'Number motif':np.array([]),
                         'p_value':np.array([]), 'fg_matches':np.array([]),
                         'fg_size':np.array([]), 'bg_size':np.array([]),
                         'bg_matches':np.array([])}) 
    return result

def motif_counter(args, vector, motif, result, dataset_info, previous_info):
    acid_location, acid_number = acids_n_l(motif)
    observed, expected, fasta, p_value=chi2_motifs(args, vector, acid_location, acid_number , dataset_info, previous_info)
    if p_value!=None:
        motif_l=letter_motif(args, acid_location, acid_number,acids=utils.ACIDS_LIST)

        result=result.append({'motif':motif_l,'Number motif':motif,
                         'p_value':p_value, 'fg_matches':observed,
                         'fg_size':previous_info[0], 'bg_size':previous_info[1],
                         'bg_matches':fasta}, ignore_index = True)
    return result

   
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
            
                prev_fg_occ = previous_table['fg_matches'][i]
                prev_bg_occ = previous_table['bg_matches'][i]
                previous_info = prev_fg_occ, prev_bg_occ
              
                result = motif_counter(args,n_vector,motif,result,dataset_info,previous_info)
                        
    b=len(result['fg_matches'].values)
    table = (result.loc[result['p_value'] < args.p_value / b][result['fg_matches']>=args.occurrences]).reset_index()
    return table           
    


def get_motifs(args, vector, previous_motifs_table, dataset_info, step):
    previous_vector = n_vectors(args, previous_motifs_table)
    matrix = multiplying_matrix(args, previous_vector, vector)
    result = table_creator()
    result = n_motifs_result(args, previous_vector, vector, previous_motifs_table, matrix, result, dataset_info)
    result.drop_duplicates(subset = ['motif'], inplace=True)
    logging.info("%d - AA motifs:\n%s" % (step, result['motif']))

    return result

def chi2_alg(fg_intervals, bg_intervals, args, results_saving_dir):
    fg_occ = utils.get_occurences(fg_intervals)
    bg_occ = utils.get_occurences(bg_intervals)
    p_values = calculate_p_value(fg_occ, bg_occ, len(fg_intervals), len(bg_intervals), args.interval_length)
    vector, single = primary_motifs(fg_occ, bg_occ ,p_values, args)
    double = double_motifs(vector, (fg_intervals, bg_intervals), args, single)   
    result_table = pd.concat([single, double], ignore_index=True)
    previous_motifs_table = double
    step = 3
    while True:
        new_motifs_table = get_motifs(args, vector, previous_motifs_table, (fg_intervals, bg_intervals), step)

        if new_motifs_table.empty:
            break
        else:                  
            result_table = pd.concat([result_table, new_motifs_table], ignore_index=True)
            previous_motifs_table = new_motifs_table
            step += 1
    cols = ['motif', 'p_value', 'fg_matches', 'fg_size', 'bg_size', 'bg_matches']     
    result_table.loc[:, cols].to_csv(os.path.join(results_saving_dir, 'motifs.csv'),index=False)            
