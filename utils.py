# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:37:30 2020

@author: Lab006

"""

ACIDS_LIST = ['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D','-']
import os
from urllib.request import urlretrieve
import pylab
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyteomics import fasta
import gzip
#import seaborn as sns
# from scipy.stats import binom
import math
import random
import profile
from scipy.stats import chisquare

import argparse
import logging
import chi2
import binomial
from collections import defaultdict, Counter
import re

def saving(args):

    sample_saving_dir = os.path.join(args.working_dir, args.name_sample)
    results_saving_dir = os.path.join(sample_saving_dir, args.modification + '_' +
                                      args.modification_site + str(args.interval_length) + '_' + args.algorithm)

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
# def peptide_identification(path_FASTA,Peptides,sample_saving_dir):
# #    print ('Identification of peptides')
#     path_identificated_peptides=os.path.join(sample_saving_dir,'peptide_identification.csv')
#     if not os.path.isfile(path_identificated_peptides):
#         FASTA_dict=dict()
#         for description, sequence in fasta.read(open(path_FASTA)):
#             FASTA_dict[sequence]=description
#         for elem in Peptides.index:
# #            print(elem)
#             #выделили элемент, будем искать его в FASTA
#             for protein in FASTA_dict:
#                 if (protein.count(elem)!=0) and (Peptides['Protein'][elem] is not None):
#                     Peptides['Protein'][elem]=None
#                 if (protein.count(elem)!=0) and (Peptides['Protein'][elem] is None):
#                     Peptides['Protein'][elem]=protein            
#         idPeptides=Peptides.dropna()
        
#         idPeptides.to_csv(path_identificated_peptides,mode='w')
#     else:
#         idPeptides=pd.read_csv(path_identificated_peptides, sep=',')
        
# #    print ('Peptides are identificated')
#     logging.info(msg=str(len(idPeptides['Protein'].values))+' of '+str(len(Peptides['Protein'].values))
#     +' modified peptides are identificated')
#     logging.debug(msg=u'Peptides are identificated')                                
#     return idPeptides


#нарезка интервалов
# def interval_maker_experimental(row, bg_fasta, args):
#     pass
    # for re.finditer('*', row['Mod_Peptide']):
    # intervals=[]
    #фиксируем пептид с модификациями
    #print('!',idPeptides)
    # for elem in idPeptides.index:
    #     peptide,mod_peptide,protein = idPeptides['Peptide'][elem],idPeptides['Mod_Peptide'][elem],idPeptides['Protein'][elem]

    #     #число модификации в пептиде
    #     asterisk=0
    #     #бежим по выбранному пептиду
    #     for k in range(len(mod_peptide)):
    #         #проверяем, есть ли модификация
    #         if mod_peptide[k]=='*':

    #             pep_site=k-1-asterisk #позиция модификации в пептиде

    #             asterisk+=1


    #             #для каждой модификации нужно выделить интервал
    #             left=pep_site-interval_length
    #             right=pep_site+interval_length


    #             #возможно,этот интервал можно вырезать из нашего пептида, проверим эту гипотезу
    #             if ((left-interval_length)>=interval_length) and (right<=len(peptide)-1):
    #                 interval=peptide[left:right+1]
    #                 intervals.append(interval)
    #             #если так не работает, то нужно обратиться к белку
    #             else:
    #                 #добавим проверку на наличие символа '\n' в конце белка
    #                 if protein[-1]=='\n':
    #                     protein=protein[:-1]
    #                 a=protein.find(peptide)
    #                 protein_site=a+pep_site #позиция модификации в белке

    #                 left=protein_site-interval_length
    #                 right=protein_site+interval_length


    #                 #возможно,этот интервал можно вырезать из нашего белка, проверим эту гипотезу
    #                 if (left>=0) and (right<=len(protein)-1):
    #                     interval=protein[left:right+1]
    #                     intervals.append(interval)
                    
    #                 #рассматриваем случай, когда сайт прибит к N-концу белка
    #                 elif ((left<0) and (right<=len(protein)-1)):
    #                     interval=' '*(interval_length-protein_site)+protein[:right+1]
    #                     intervals.append(interval)
                        
    #                 #рассматриваем случай, когда сайт прибит к C-концу белка        
    #                 elif ((left>=0) and (right>len(protein)-1)):
    #                     interval=protein[left:]+' '*(right+1-len(protein))
    #                     intervals.append(interval)
                    
    #                 #рассматриваем случай, когда белок маленький 
    #                 else:
    #                     interval=' '*(interval_length-protein_site)+protein+' '*(right+1-len(protein))
    #                     intervals.append(interval)

    # mod_intervals=[]

    # for interval in intervals:
    #     if interval[interval_length]==modification_site:
    #         mod_intervals.append(interval)                   
    # return mod_intervals

def fasta_match(row, bg_fasta, interval_length,modification_site):
    intervals = []
    k=0
    print(row['Peptide'],k)
    for name, seq in bg_fasta.items():
        i = 0
        start = seq[i:].find(row['Peptide'].replace('*',''))
        i = start
        while start >= 0:
            for asterisks, modif in enumerate(re.finditer('\*', row['Peptide']), 1):
                interval_start = i + modif.span()[0] - interval_length - asterisks
                interval_end = interval_start + 2 * interval_length + 1
                interval=seq[interval_start: interval_end]
                if interval[interval_length]==modification_site:
                    intervals.append(interval)
                else:
                    print('wrong',interval)
            start = seq[i+1:].find(row['Peptide'].replace('*',''))
            i += start + 1
    k+=1        
    return intervals

#формирование набора интервалов
# def mod_intervals_DB_experimental(Peptides, args, bg_fasta, sample_saving_dir):
# #    print('Making intervals DB')
#     # idPeptides=peptide_identification(args.fasta, Peptides, sample_saving_dir)


#     mod_intervals=interval_maker_experimental(idPeptides, args.interval_length, args.modification_site)
# #    print('Intervals DB is ready!')
#     logging.debug(u'Intervals DB is ready')
#     logging.info(u'Set of {} modified intervals is ready'.format (len(mod_intervals))) 
#     return mod_intervals

# #правильная функция
# def FASTA_DB_creation(path_FASTA):
#     background_DB=dict()
#     #берем белок из белковой БД человека
#     for description, sequence in fasta.read(open(path_FASTA)):    
#         #ищем его номер
#         i=description.find('|')
#         k=description.rfind('|')
#         s=description[i+1:k]
#         #если белка нет в background БД записываем 
#         if s not in background_DB:
#             if sequence.count('X')==0:
#                 background_DB[s]=sequence           
#     return  background_DB    


# #функция выделения интервалов для составления background
# def Background_creator(elem,i,interval_length):
#     #определим интервал поиска аминокислот        
#     #рассматриваем случай, когда сайт прибит к N-концу белка
#     if (i<=interval_length-1) and (len(elem)-i>interval_length):
#         #количество кислот слева
#         left=i
#         #количество рассматриваемых кислот справа
#         right=interval_length
#         interval=' '*(interval_length-left)+elem[i-left:i+right+1]
#     #рассматриваем случай, когда сайт прибит к C-концу белка    
#     elif  (i>interval_length-1) and (len(elem)-i<=interval_length):
#         left=interval_length
#         right=len(elem)-i-1
#         interval=elem[i-left:i+right+1]+' '*(interval_length-right)
#     #рассматриваем случай, когда белок маленький    
#     elif (i<=interval_length-1) and (len(elem)-i<=interval_length):
#         left=i
#         right=len(elem)-i-1
#         interval=' '*(interval_length-left)+elem[i-left:i+right+1]+' '*(interval_length-right)
#     #случай, когда можем рассмотреть интервал (-15,+15)    
#     else:
#         right=left=interval_length
#         interval=elem[i-left:i+right+1]
#     #выделяем наш интервал-последовательность из белка
#     #interval=elem[number-left:number+right+1]
#     return interval


def background_maker(args):
#    print('Making background DB')
    #хотим сделать background из идентифицированных белков
    bg_fasta = dict()
    bg = defaultdict()
    background = set()
    with fasta.read(args.fasta) as f:
        for name, sequence in f:
            name_id = name.split('|')[1]
            extended_seq = ''.join(['-' * args.interval_length, sequence, '-' * args.interval_length])
            bg_fasta[name_id] = extended_seq
            mod_aa_indexes = re.finditer(args.modification_site, extended_seq)
            bg_intervals = [extended_seq[i.span()[0] - args.interval_length: i.span()[0] + args.interval_length + 1] for i in mod_aa_indexes]
            bg[name_id] = bg_intervals
            background.update(bg_intervals)

    logging.info(u'Set of ' + str(len(background))+ u' background intervals is created')
    logging.debug(u'Background DB is ready')    
    return background, bg_fasta   


#функции для валидации


#N_calculator= lambda intervals:len(intervals)

def aa_counter(col):
    return pd.Series(Counter(col))

def get_occurences(intervals_list, interval_length, saving_file, acids=ACIDS_LIST):
    #создаем матрицу для подсчета числа n и забиваем ее нулями
    # n=[]
    # for i in range(len(acids)):
    #     a=[0]*(interval_length*2+1)
    #     n.append(a)
    # for interval in intervals:
    #     for k in range(len(interval)):
    #         #проверяем какая кислота из заданного набора стоит на зафиксированной позиции
    #         for i in range(len(acids)):
    #             #номер аминокислоты i в списке аминокислот совпадает с номером этой аминокислоты в матрице n
    #             acid=acids[i]
    #             #если совпадение есть,то добавляем единичку в матрицу n на соответствующую позицию
    #             if interval[k]==acid:
    #                 n[i][k]+=1
    #             else:
    #                 continue
    # path=os.path.join(results_saving_dir,'occurrences.txt')
    # saving=open(path,'w')
    # for i in range(len(acids)):
    #     for k in range(interval_length*2+1):
    #         saving.write(str(n[i][k])+' ')
    #     saving.write('\n')
    # saving.close()
    # logging.debug(msg=u'Occurrence matrix was created')                    
    # return n
    # occurance = pd.DataFrame(index=acids)
    # i=((intervals_list.replace('[','')).replace("'",'')).split(']')
    # print(i[:-1])
    df = pd.DataFrame([list(i) for i in intervals_list], columns=range(-interval_length, interval_length + 1))
    occ = df.apply(aa_counter, axis=0)
    print(occ)
    occ.to_csv(saving_file, sep='\t')
    return occ

# def background_n_matrix(interval_length,background,results_saving_dir, acids=ACIDS_LIST):
#     background_n=[]
#     path=os.path.join(results_saving_dir,'background.txt')
#     saving=open(path,'w')
#     #забиваем матрицу нулями
#     for i in range(len(acids)):
#         a=[0]*(interval_length*2+1)
#         background_n.append(a)
#     #на данном этапе считаем встречаемость аминокислоты на j позиции
#     for interval in background:
#     #фиксируем позицию в интервале
#         for k in range(len(interval)):
#             #проверяем какая кислота из заданного набора стоит на зафиксированной позиции
#             for i in range(len(acids)):
#                 #номер аминокислоты i в списке аминокислот совпадает с номером этой аминокислоты в матрице n
#                 acid=acids[i]
#                 #если совпадение есть,то добавляем единичку в матрицу n на соответствующую позицию
#                 if interval[k]==acid:
#                     background_n[i][k]+=1
#                 else:
#                     continue
              
    # for i in range(len(acids)):
    #     for k in range(interval_length*2+1):
    #         if k!=interval_length:
    #             saving.write(str(background_n[i][k])+' ')
    #         else:
    #             saving.write(str(0)+' ')
    #     saving.write('\n')
    # saving.close()
    # logging.debug(msg=u'Background occurrence matrix was created')                     
    # return background_n

def saving_table(results_saving_dir,result,interval_length,name):
    path=os.path.join(results_saving_dir,'table'+str(interval_length)+'_'+name+'.csv')
    result.to_csv(path)   
    
def peptides_table(args,sample_saving_dir,bg_fasta):
    if os.path.exists(os.path.join(sample_saving_dir, 'peptide_identification.csv')):
        
        idPeptides=pd.read_csv(os.path.join(sample_saving_dir, 'peptide_identification.csv'),sep=',',usecols=[1,2])
        for i in range(len(idPeptides['fasta_match'].index)):
            idPeptides['fasta_match'][i]=(((idPeptides['fasta_match'][i].replace('[','')).replace(']','')).replace("'",'')).split(',')
    else:    
        peptides=[]
        with open(args.dataset) as f:
            for line in f:
                peptides.append(line[:-1])
        Peptides=pd.DataFrame({'Peptide':peptides})
        Peptides['fasta_match'] = Peptides.apply(fasta_match, args=[bg_fasta, args.interval_length, args.modification_site], axis=1)
        # print(background)np
        Peptides['unique'] = Peptides.apply(lambda x: True if len(x['fasta_match']) == 1 else False, axis=1)
        # print(Peptides)
        idPeptides = Peptides[Peptides['unique'] == True]
        idPeptides.to_csv(os.path.join(sample_saving_dir, 'peptide_identification.csv'), mode='w')
        
    return idPeptides    

def output(args):
    sample_saving_dir,results_saving_dir = saving(args)
    background, bg_fasta = background_maker(args)
    idPeptides=peptides_table(args,sample_saving_dir,bg_fasta)
    print(idPeptides)
    # Peptides['fasta_match'] = Peptides.apply(fasta_match, args=[bg_fasta, args.interval_length], axis=1)
    # # print(background)
    # Peptides['unique'] = Peptides.apply(lambda x: True if len(x['fasta_match']) == 1 else False, axis=1)
    # # print(Peptides)
    # idPeptides = Peptides[Peptides['unique'] == True]
    # print(idPeptides)
    
    # idPeptides.to_csv(os.path.join(sample_saving_dir, 'peptide_identification.csv'), mode='w')

    # intervals = mod_intervals_DB_experimental(Peptides, args, bg_fasta, sample_saving_dir)
    # intervals = interval_maker_experimental(idPeptides, args, bg_fasta)
    if args.algorithm=="binom":
#        print('I AM HERE')
        P = binomial.P_counter_bi(intervals, args.interval_length, args.modification_site, background, results_saving_dir, acids=ACIDS_LIST)
        occurrences = binomial.occurrences_counter_bi(intervals, args.interval_length, args.modification_site, results_saving_dir, acids=ACIDS_LIST)
        P_final = binomial.final_validation_bi(args.interval_length, occurrences, P, acids=ACIDS_LIST)
        single, double, triple, quadruple = binomial.motifs_bi(P_final, args.interval_length, args.modification_site, 
                                                               args.background, intervals, P, results_saving_dir, acids=ACIDS_LIST)

        logging.info(msg='Program was finished successfully') 
        return P, occurrences, intervals, background
    else:
        
        occurrences = get_occurences( (idPeptides['fasta_match']).sum(), args.interval_length,
                                     os.path.join(results_saving_dir, 'occurences.csv'), 
                                     acids=ACIDS_LIST)
        background_n = get_occurences(background, args.interval_length, 'background.csv')
#        expected_FASTA, chi2_results, chi2_selection = chi2.p_value(background_n, occurrences, args.interval_length,
#                                                                    args.modification_site,results_saving_dir, acids=ACIDS_LIST)
        p_value=chi2.p_value(occurrences,background_n,args.interval_length,results_saving_dir)
#        single, double, triple, quadruple = chi2.motifs(chi2_selection, args.interval_length, args.modification_site, background,
#                                                        intervals, results_saving_dir, acids=ACIDS_LIST)
#        single, double, triple, quadruple = chi2.motifs(p_value, args.interval_length, args.modification_site, background_n,
#                                                occurrences, results_saving_dir, acids=ACIDS_LIST)
        vector,table=chi2.motifs(occurrences,background_n,p_value,args,results_saving_dir)
        logging.info(msg='Program was finished successfully') 
#        return chi2_results,chi2_selection,intervals,background
        return vector,table

