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
import chi2
import binom

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
    logging.info(u'Set of '+str(len(mod_intervals))+u' modified intervals is ready') 
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
    expected_FASTA,chi2_results,chi2_selection=chi2.p_value(background_n,occurrences,interval_length,modification_site,acids,results_saving_dir)
    single,double,triple,quadruple=chi2.motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,results_saving_dir)
#    print(single,double,triple,quadruple)
#    logging.info(msg='Program was finished successfully')

    
    return chi2_results,chi2_selection,intervals,background
    
def output_experimental_bi(working_dir,name_sample,Peptides,interval_length,modification_site,modification,path_FASTA,algorithm):
    sample_saving_dir,results_saving_dir=saving(working_dir,name_sample,interval_length,modification,modification_site,algorithm)
    acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
    intervals=mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,sample_saving_dir)
#    print('Длина интервалов в exp',len(intervals))
    background=background_maker(modification_site,interval_length,path_FASTA,modification)
#    print('Длина интервалов в back',len(background))
    P=binom.P_counter_bi(intervals,interval_length,modification_site,acids,background,results_saving_dir)
    occurrences=binom.occurrences_counter_bi(intervals,interval_length,acids,modification_site,results_saving_dir)
    P_final=binom.final_validation_bi(acids,interval_length,occurrences,P)
    single,double,triple,quadruple=binom.motifs_bi(acids,P_final,interval_length,modification_site,background,intervals,P,results_saving_dir)
#    print(single,double,triple,quadruple)
#    logging.info(msg='Program was finished successfully')
    return P,occurrences,intervals,background    

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
