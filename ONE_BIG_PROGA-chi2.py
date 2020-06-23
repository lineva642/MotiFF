
import os
from urllib.request import urlretrieve
import pylab
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyteomics import fasta
import gzip
import seaborn as sns
from scipy.stats import binom
import math
import random
import profile
from scipy.stats import chisquare



#набор функции для построения базы модифицированных интервалов
#для экспериментальных наборов данных



#оптимизация функции peptide_identification(path_FASTA,Peptides,path_sample)
def peptide_identification(path_FASTA,Peptides,path_sample):
    print ('Identification of peptides')
    path_identificated_peptides=path_sample+'/'+'peptide_identification.csv'
    if not os.path.isfile(path_identificated_peptides):
        FASTA_dict=dict()
        for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):
            FASTA_dict[sequence]=description
        for elem in Peptides.index:
            print(elem)
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
    print ('Peptides are identificated')
    return idPeptides


#для Иры
def peptide_identification(path_FASTA,Peptides,path_sample):
    print ('Identification of peptides')
    path_identificated_peptides=path_sample+'/'+'peptide_identification.csv'
    if not os.path.isfile(path_identificated_peptides):
        FASTA_dict=dict()
        for description, sequence in fasta.read('/home/vikalin/Ira/sp_human_canonical_21032019_reversed.fasta'):
            FASTA_dict[sequence]=description
        for elem in Peptides.index:
            print(elem)
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
    print ('Peptides are identificated')
    return idPeptides


#улучшение программы
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


#улучшение проги
def mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,path_sample):
    print('Making intervals DB')
    #proteins,peptides,mod_peptides=peptide_identification(path_FASTA,Peptides,path_sample)
    idPeptides=peptide_identification(path_FASTA,Peptides,path_sample)
    
    #mod_intervals=interval_maker_experimental(idPeptides['Protein'].values,idPeptides['Peptide'].values,idPeptides['Mod_Peptide'].values,interval_length,modification_site)
    mod_intervals=interval_maker_experimental(idPeptides,interval_length,modification_site)
    print('Intervals DB is ready!')
    return mod_intervals


# In[16]:


#второй датасет
def mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,path_sample):    
    #mod_intervals=interval_maker_experimental(idPeptides['Protein'].values,idPeptides['Peptide'].values,idPeptides['Mod_Peptide'].values,interval_length,modification_site)
    mod_intervals=interval_maker_experimental(Peptides,interval_length,modification_site)
    print('Intervals DB is ready!')
    return mod_intervals


# In[48]:


#набор функции для построения базы модифицированных интервалов
#в случае, когда у нас в качестве входных данных рассматривается ПТМ


# In[90]:


def PTM_DB_initial(path_PTM,modification,organism):
    #в случае, когда у нас в качестве входных данных рассматривается ПТМ
    #path_PTM='/home/vikalin/knowledgebase/PTM/PTM.txt'
    #modification='ACETYLATION'
    #organism='Human'

    initialDB=open(path_PTM,'r')
    names=[]
    sites=[]
    k=0
    for line in initialDB:
        if line.count(modification)!=0 and line.count(organism)!=0:
            print(k)
            k+=1
            l=line.split('\t')
            #выделяем из этого массива номер белка и сайт модификации
            names.append(l[2])
            sites.append(l[5])
    Proteins_sites=pd.DataFrame({'Sites':sites}, index=names)
    Proteins_sites['Proteins']=None
    initialDB.close()
    return Proteins_sites


# In[91]:


def mod_proteins_DB(path_FASTA,path_sample,Proteins_sites):
    i=0
    #на данном этапе происходит идентификация белков в базе ПТМ из FASTA
    #file=open('/home/vikalin/Article/identification_example.txt','r')
    #names=[]
    #for line in file:
    #    l=line.split('\t')
    #    names.append(l[0])
    #file.close()
    #file=open('/home/vikalin/Article/identification_example.txt','a')
    path_identificated_proteins=path_sample+'/'+'protein_identification.txt'
    if not os.path.isfile(path_identificated_proteins):
        identification=open(path_identificated_proteins,'w')
        for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):
            s=description[description.find('|')+1:description.rfind('|')]
            if (s in Proteins_sites.index):
                print(i)
                i+=1
                Proteins_sites.loc[s,'Proteins']=sequence
                identification.write(s+'/t'+sequence+'\n')
            else:
                print('NO',i)
        print('THE END')
        
    else:
        identification=open(path_identificated_proteins,'r')
        for line in identification:
            l=line.split('/t')
            number=l[0]
            protein=l[1][:-1]
            if (number in Proteins_sites.index):
                print(i)
                i+=1
                Proteins_sites.loc[number,'Proteins']=protein
            else:
                print('NO',i)
        print('THE END')    
        
    Proteins_sites=Proteins_sites.dropna()
    Proteins_sites.index=[i for i in range(len(Proteins_sites))]
    mod_Proteins=Proteins_sites.copy()
    print('THE END 1')
    l=0
    for i in range(len(mod_Proteins)):
        protein=mod_Proteins.loc[i,'Proteins']
        site=mod_Proteins.loc[i,'Sites']
        number=int(site[1:])-1
        acid=site[0]
        if number<len(protein):
            if protein[number]==acid:
                print(l)
                l+=1
                mod_protein=protein[:number+1]+'*'+protein[number+1:]
                mod_Proteins.loc[i,'Proteins']= mod_protein
        else:
            print(protein,number)
    del mod_Proteins['Sites']
    identification.close()
    print('THE END 2')
    mod_Proteins.to_csv(path_sample+'/'+'identification.csv')
    return mod_Proteins 


# In[92]:


def interval_maker_FASTA(mod_Proteins,interval_length,modification_site):
    print('YES!')
    intervals=[]
    proteins=mod_Proteins['Proteins'].values
    for protein in proteins:
        site=protein.find('*')-1 #позиция модификации в белке
        if protein[site]==modification_site:
            protein=protein.replace('*','')
            #для каждой модификации нужно выделить интервал
            left=site-interval_length
            right=site+interval_length
           # print(protein,left,right)
            #возможно,этот интервал можно вырезать из нашего белка полнотью, проверим эту гипотезу
            if ((left>=0) and (right<=len(protein)-1)):
                interval=protein[left:right+1]
                intervals.append(interval)
            #рассматриваем случай, когда сайт прибит к N-концу белка
            elif ((left<0) and (right<=len(protein)-1)):
                interval=' '*(interval_length-site)+protein[:right+1]
                intervals.append(interval)

            #рассматриваем случай, когда сайт прибит к C-концу белка        
            elif ((left>=0) and (right>len(protein)-1)):
                interval=protein[left:]+' '*(right+1-len(protein))
                intervals.append(interval)

            #рассматриваем случай, когда белок маленький 
            else:
                interval=' '*(interval_length-site)+protein+' '*(right+1-len(protein))
                intervals.append(interval)
                 
    return intervals


# In[93]:


def mod_intervals_DB_PTM(path_PTM,path_FASTA,path_sample,modification,organism,interval_length,modification_site):
    #случай, когда идентификации не было, таблица модифицированных белков не составлена
    if not os.path.isfile(path_sample+'/'+'identification.csv'):
        Proteins_sites=PTM_DB_initial(path_PTM,modification,organism)
        mod_Proteins=mod_proteins_DB(path_FASTA,path_sample,Proteins_sites)
    #случай, когда идентификация была, таблица модифицированных белков составлена    
    else:
        mod_Proteins = pd.read_csv(path_sample+'/'+'identification.csv', sep=',')
    intervals=interval_maker_FASTA(mod_Proteins,interval_length,modification_site)
    print('Количество интервалов из PTM базы:', len(intervals))
    return intervals,len(intervals)


# In[94]:


#набор функции для построения базы немодифицированных интервалов или backgroundов


# In[95]:


def PTM_DB_creation(path_PTM,path_FASTA,modification):
    #path_PTM-'/home/vikalin/knowledgebase/PTM/PTM.txt'
    initialDB=open(path_PTM,'r')
    path_PTM_mod='/home/vikalin/knowledgebase/PTM/human_'+modification+'.txt'
    humanDB=open(path_PTM_mod,'w')
    for line in initialDB:
        if line.count(modification)!=0 and line.count('Human')!=0:
            humanDB.write(line)
    initialDB.close()
    humanDB.close()        
    humanDB=open(path_PTM_mod,'r')
    PTM=dict()
    for line in humanDB:
        #представляем строку из БД в виде массива
        l=line.split('\t')
        #выделяем из этого массива номер белка и сайт модификации
        name=l[2]
        site=l[5]

        if name in PTM:
            PTM[name].add(site)
        else:
            PTM[name]=set()
            PTM[name].add(site)
    humanDB.close()  
    #создаем FASTA базу,где description-номер белка и аминокислотная последовательность, sequence-последовательность сайтов
    #PHOSPHORYLATION_HUMAN=open('/home/vikalin/Article/PHOSPHORYLATION_HUMAN_DB.fasta','w')
    PTM_DB=dict()
    #берем белок из белковой БД человека
    #path_FASTA-'/home/vikalin/knowledgebase/PTM/HUMAN_proteinDB.fasta.gz'
    for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):

        #ищем его номер
        i=description.find('|')
        k=description.rfind('|')
        s=description[i+1:k]
        #если белок есть в словаре PTM,то должны записать его в новый словарь PHOSPHORYLATION, где его аминокислотная последовательность-ключ,сайты модификации-значение
        if s in PTM:
            #запишем его множество сайтов в строку
            sites=''
            for site in PTM[s]:
                if sites=='':
                    sites=str(site)
                else:    
                    sites=sites+' '+str(site)       
            #запишем эту строку в новый словарь

            #нужно исключить белки с аминокислотой 'X'
            if (sites!='') and (sequence.count('X')==0):
                PTM_DB[sequence]=sites
    #            PHOSPHORYLATION_HUMAN.write('>'+sequence+'\n'+sites+'\n')
        else:
            continue
    #PHOSPHORYLATION_HUMAN.close()
    return PTM_DB


# In[13]:


#правильная функция
def FASTA_DB_creation(path_FASTA):
    background_DB=dict()
    #берем белок из белковой БД человека
    for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):    
        #ищем его номер
        i=description.find('|')
        k=description.rfind('|')
        s=description[i+1:k]
        #если белка нет в background БД записываем 
        if s not in background_DB:
            if sequence.count('X')==0:
                background_DB[s]=sequence           
    return  background_DB    


# In[55]:


#для рандома
def FASTA_DB_creation(path_FASTA):
    background_DB=dict()
    #берем белок из белковой БД человека
    for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):    
        #ищем его номер
        i=description.find('|')
        k=description.rfind('|')
        s=description[i+1:k]
        #если белка нет в background БД записываем 
        if s not in background_DB:
            if sequence.count('X')==0:
                str_var = list(sequence)
                np.random.shuffle(str_var)
                string=''.join(str_var)
                background_DB[s]=string           
    return  background_DB 


# In[96]:


#для HeLa
def FASTA_DB_creation(path_FASTA):
    proteins=pd.read_csv('/home/vikalin/20151001_03_Qp1_HeLa_1ug.mzML.demix_proteins.csv',sep='\t')
    hela_protein=proteins['dbname'].values
    hela=dict()
    for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):    
        #ищем его номер
        i=description.find('|')
        k=description.rfind('|')
        s=description[i+1:k]
        if (s in hela_protein) and sequence.count('X')==0:
            hela[s]=sequence
    return  hela       


# In[13]:


#для Иры
def FASTA_DB_creation(path_FASTA):
    background_DB=dict()
    #берем белок из белковой БД человека
    for description, sequence in fasta.read('/home/vikalin/Ira/sp_human_canonical_21032019_reversed.fasta'):    
        #ищем его номер
        i=description.find('|')
        k=description.rfind('|')
        s=description[i+1:k]
        #если белка нет в background БД записываем 
        if s not in background_DB:
            if sequence.count('X')==0:
                background_DB[s]=sequence           
    return  background_DB    


# In[97]:


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


# In[98]:


def background_array_PTM(modification_site,interval_length,path_PTM,path_FASTA,modification):
    #хотим сделать background из идентифицированных белков
    PTM_DB=PTM_DB_creation(path_PTM,path_FASTA,modification)
    background=[]
    for elem in PTM_DB:
        #выделяем строку сайтов модификации
        s=PTM_DB[elem]
        #сделаем массив для хранения сайтов в виде объектов
        sites=s.split(' ')
        #бежим по аминокислотной последовательности
        for i in range(len(elem)):
            #ищем совпадающие с сайтом модификации
            if elem[i]==modification_site:
                #восстанавливаем сайты модификации
                number=elem[i]+str(i+1)
                #проверяем, является ли аминокислота модифицированной
                if number in sites:
                    #exaption.append(number)
                    continue
                else:
                    interval=Background_creator(elem,i,interval_length)
                    background.append(interval)
    return background 


# In[99]:


def background_array_FASTA(path_FASTA,modification_site,interval_length,intervals):
    #хотим сделать background из идентифицированных белков
    background_DB=FASTA_DB_creation(path_FASTA)
    print('background_array_FASTA',len(background_DB))
    background=[]

    for elem in background_DB:
        s=background_DB[elem]
        for i in range(len(s)):
            #ищем совпадающие с сайтом модификации
            if s[i]==modification_site:
                interval=Background_creator(s,i,interval_length)
                background.append(interval)

          
#    random.shuffle(background)
#    print('I shuffled background!')
   # print('Previous background intervals:',len(background))
#    background=background[:len(intervals)]


   # print('New background intervals:',len(background))           
    return background 


# In[100]:


def background_maker(modification_site,interval_length,path_FASTA,modification,intervals,path_PTM=None):
    print('path PTM',path_PTM)
    print('Making background DB')
    if path_PTM==None:
        print('I am doing right thing!')
        background=background_array_FASTA(path_FASTA,modification_site,interval_length,intervals)
    else:
        print('I am doing bad thing!')
        background=background_array_PTM(modification_site,interval_length,path_PTM,path_FASTA,modification)
    print('Background DB is ready')    
    return background   


# In[101]:


#функции для валидации


# In[102]:


N_calculator= lambda intervals:len(intervals)


# In[103]:


def n_calculator(acids,intervals,interval_length,path_results):
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
    path=path_results+'/'+'occurrences'+'.txt'
    saving=open(path,'w')
    for i in range(len(acids)):
        for k in range(interval_length*2+1):
            saving.write(str(n[i][k])+' ')
        saving.write('\n')
    saving.close()                    
    return n


# In[104]:


def background_n_matrix(acids,interval_length,background,path_results):
    background_n=[]
    path=path_results+'/'+'background'+'.txt'
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
    return background_n


# In[105]:


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
    return expected_FASTA,all_exp
#    return background_n,all_exp


# In[106]:


# считаем p-value
def chi2_result(occurrences_FASTA,expected_FASTA,all_exp,interval_length,path_results):
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
    path=path_results+'/'+'p_value'+'.txt'
    saving=open(path,'w')
    for i in range(0,21):
        for k in range(interval_length*2+1):
            saving.write(str(chi2_results[i][k])+' ')
        saving.write('\n')
    saving.close()
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
    return chi2_selection            


# In[108]:


#напишем одну функцию для подсчета p-value


# In[109]:


def p_value(background_n,occurrences,interval_length,modification_site,acids,path_results):
    expected_FASTA,all_exp=expected_destribution(background_n,occurrences,interval_length)
    chi2_results=chi2_result(occurrences,expected_FASTA,all_exp,interval_length,path_results)
    chi2_selection=p_value_selection(chi2_results,occurrences,interval_length)
    heatmap_visualization(chi2_selection,acids,interval_length,modification_site,'Отбор по p-value',path_results,name='p_value_selection')
    return expected_FASTA,chi2_results,chi2_selection


# In[110]:


#функции для визуализации


# In[111]:


def heatmap_visualization(matrix,acids,interval_length,modification_site,title,path_results,name=None):
    x=[]
    for i in range(-interval_length,interval_length+1):
        if i==0:
            x.append(modification_site)
        else:    
            x.append(i)
    fig = plt.figure()
    ax = sns.heatmap(matrix,xticklabels=x,yticklabels=acids,linewidths=.5,cmap="RdBu_r",center=1,square=True)
    fig.set_figwidth(15)    #  ширина и
    fig.set_figheight(15)#  высота "Figure"
    ax.set_title(title)
    
    #сохранение рисунка
    if name!=None:
        path=path_results+'/'+name+'_'+modification_site+str(interval_length)+'.png'
        plt.savefig(path)
    plt.show()


# In[112]:


#функция для подсчета комлексных мотивов


# In[113]:


def primary_motifs(acids,interval_length,chi2_selection,path_results,modification_site):
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
        print(i)
        motif=dict()
        pair_1_x,pair_1_y=primary_motif['Location'].values[i]
        motif[pair_1_y]=acids[pair_1_x]
        motif[interval_length]=modification_site.lower()
        list_keys = list(motif.keys())
        list_keys.sort()
        word_motif=motif[list_keys[0]]+'.'*(list_keys[1]-list_keys[0]-1)+motif[list_keys[1]]
        motifs.append(word_motif)
        
    primary_motif['Motifs']=motifs
    saving_table(path_results,primary_motif,interval_length,'primary')
    return primary_motif


# In[73]:


#для Иры
def primary_motifs(acids,interval_length,chi2_selection,chi2_results,path_results,modification_site):
    primary_motifs=[]
    primary_location=[]
    primary_p_value=[]
    for i in range(0,21):
        for k in range(0,2*interval_length+1):
            if chi2_selection[i][k]==1:
                primary_motifs.append(acids[i])
                primary_location.append((i,k))
                primary_p_value.append(chi2_results[i][k])
    primary_motif=pd.DataFrame({'Acid':primary_motifs,'Location':primary_location,'p-value':primary_p_value})
    
    motifs=[]
    for i in range(len(primary_motif['Location'].values)):
        print(i)
        motif=dict()
        pair_1_x,pair_1_y=primary_motif['Location'].values[i]
        motif[pair_1_y]=acids[pair_1_x]
        motif[interval_length]=modification_site.lower()
        list_keys = list(motif.keys())
        list_keys.sort()
        word_motif=motif[list_keys[0]]+'.'*(list_keys[1]-list_keys[0]-1)+motif[list_keys[1]]
        motifs.append(word_motif)
        
    primary_motif['Motifs']=motifs
    del primary_motif['Acid']
    del primary_motif['Location']
    saving_table(path_results,primary_motif,interval_length,'primary')
    return primary_motif


# In[114]:


def double_motifs(primary_motif,acids,intervals,background,path_results,interval_length,modification_site):
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
        print(i)
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
        print(i)
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
    
    saving_table(path_results,second_motifs_selection_p_copy,interval_length,'double')
    
    return second_motifs_selection_p_copy


# In[115]:


def triple_motifs(primary_motif,second_motifs,acids,intervals,background,path_results,modification_site):
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
        print(i)
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
        print(i)
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
    
    saving_table(path_results,triple_motifs_selection_p_copy_copy,interval_length,'triple')
    
    return triple_motifs_selection_p_copy_copy


# In[116]:


#делаем прогу для 4хбуквенного мотива
def quadruple_motifs(primary_motif,triple_motifs,acids,intervals,background,path_results,modification_site):
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
        print(i)
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
        print(i)
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
    
    saving_table(path_results,quadro_motifs_selection_p_copy_copy,interval_length,'quadruple')
    
    return quadro_motifs_selection_p_copy_copy


# In[117]:


#улучшенная версия 2
def motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,path_results):
    single=primary_motifs(acids,interval_length,chi2_selection,path_results,modification_site)
    double=double_motifs(single,acids,intervals,background,path_results,interval_length,modification_site)
    if double is not None:
        triple=triple_motifs(single,double,acids,intervals,background,path_results,modification_site)
        if triple is not None:
            quadruple=quadruple_motifs(single,triple,acids,intervals,background,path_results,modification_site)
        else:
            quadruple=None
    else:
        triple=None
    return single,double,triple,quadruple
#    return single


# In[62]:


#для Иры
def motifs(acids,chi2_selection,chi2_results,interval_length,modification_site,background,intervals,path_results):
    single=primary_motifs(acids,interval_length,chi2_selection,chi2_results,path_results,modification_site)
#    double=double_motifs(single,acids,intervals,background,path_results,interval_length,modification_site)
#    if double is not None:
#        triple=triple_motifs(single,double,acids,intervals,background,path_results,modification_site)
#        if triple is not None:
#            quadruple=quadruple_motifs(single,triple,acids,intervals,background,path_results,modification_site)
#        else:
#            quadruple=None
#    else:
#        triple=None
#    return single,double,triple,quadruple
    return single


# In[118]:


#организация хранения результатов
def saving(name_sample,interval_length,modification,modification_site,path_PTM=None): 
    if path_PTM!=None:    
        name_folder=modification_site+'_'+str(interval_length)+'_background:PTM'
    else:
        name_folder=modification_site+'_'+str(interval_length)+'_background:FASTA'  
    path_sample='/home/vikalin/Article/'+name_sample
    path_modification='/home/vikalin/Article/'+name_sample+'/'+modification
    print(path_modification)
    path_results=path_modification+'/'+name_folder
    print(path_results)
    if not (os.path.exists(path_modification)):
        os.mkdir(path_modification)
        if not (os.path.exists(path_results)):
            os.mkdir(path_results)
    else:
        if not (os.path.exists(path_results)):
            os.mkdir(path_results)
    return  path_sample,path_modification,path_results       


# In[119]:


#улучшение
#хранение таблицы мотивов
def saving_table(path_result,result,interval_length,name):
    path=path_result+'/table'+str(interval_length)+'_'+name+'.csv'
    result.to_csv(path)   


# In[120]:


def output_experimental(name_sample,Peptides,interval_length,modification_site,modification,path_FASTA,path_PTM=None):
    #path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
    #path_identification='/home/vikalin/Article/identification.txt'
    #path_PTM='/home/vikalin/knowledgebase/PTM/human_phosphorylation.txt'
    path_sample,path_modification,path_results=saving(name_sample,interval_length,modification,modification_site)
    acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
    intervals=mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,path_modification)
    print('Длина интервалов в exp',len(intervals))
    background=background_maker(modification_site,interval_length,path_FASTA,modification,intervals,path_PTM)
    print('Длина интервалов в back',len(background))
    occurrences=n_calculator(acids,intervals,interval_length,path_results)
    background_n=background_n_matrix(acids,interval_length,background,path_results)
    expected_FASTA,chi2_results,chi2_selection=p_value(background_n,occurrences,interval_length,modification_site,acids,path_results)
    single,double,triple,quadruple=motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,path_results)
    print(single,double,triple,quadruple)
#     single,double,triple,quadruple=motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,path_results)
#     print(single,double,triple,quadruple)
    
    return chi2_results,chi2_selection,intervals,background


# In[121]:


def output_PTM(name_sample,organism,interval_length,modification_site,modification,path_FASTA,path_PTM):
    #path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
    #path_identification='/home/vikalin/Article/identification.txt'
    #path_PTM='/home/vikalin/knowledgebase/PTM/human_phosphorylation.txt'
    path_sample,path_modification,path_results=saving(name_sample,interval_length,modification,modification_site,path_PTM)
    acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
    intervals,length_intervals=mod_intervals_DB_PTM(path_PTM,path_FASTA,path_modification,modification,organism,interval_length,modification_site)
    background=background_maker(modification_site,interval_length,path_FASTA,modification,length_intervals)
    print(background[:5])
    P=P_counter(intervals,interval_length,modification_site,acids,background,path_results)
    heatmap_visualization(P,acids,interval_length,modification_site,'Исходные вероятности для '+modification+' '+modification_site,path_results,'initial_probability')
    occurrences=occurrences_counter(intervals,interval_length,acids,modification_site,path_results)
    heatmap_visualization(occurrences,acids,interval_length,modification_site,'Occurrences для '+modification+' '+modification_site,path_results,'initial_occurrences')
    volcano_plot(acids,interval_length,occurrences,P,modification_site,modification,path_results,'plot')
    P_final=final_validation(acids,interval_length,occurrences,P)
    heatmap_visualization(P_final,acids,interval_length,modification_site,'Отбор по вероятностям и встречаемости '+modification+' '+modification_site,path_results,'final_probability')
    single,double,triple=motifs(acids,P_final,interval_length,modification_site,background,intervals,P,path_results)
    print(single,double,triple)
    results=result(single,double,triple,background_type,interval_length)
    print(results)
    saving_table(path_sample,results,interval_length)
    print(results)
    return results


# In[122]:


#поверка на наборе данных


# In[123]:


#загружаем таблицу модифицированных пептидов
if not os.path.isfile('/home/vikalin/Article/dataset.xls'):
    print ('Downloading the dataset')
    urlretrieve(
        'https://www.pnas.org/content/pnas/suppl/2004/08/06/0404720101.DC1/04720Table3.xls',
        '/home/vikalin/Article/dataset.xls')
    print ('Done!')
else:
    print('XLS file for dataset already exists!')


# In[124]:


#улучшение проги
xl = pd.ExcelFile('/home/vikalin/Article/dataset.xls')
df1 = xl.parse('Sheet1')
df1.columns=df1.loc[0]
df=df1.drop(0)
Peptides=((df.loc[df['Amb?'] == 'UNIQUE']).reset_index()).loc[:,df.columns.intersection(['Peptide']) ]
indexes=[elem.replace('*','') for elem in Peptides['Peptide']]
Peptides['Mod_Peptide']=indexes  
Peptides.columns = ['Mod_Peptide','Peptide']
Peptides['Protein']=None
ind=Peptides['Peptide'].values
Peptides['index']=ind
Peptides=Peptides.set_index('index')
Peptides


# In[144]:


#Ира первый файл
ex_1=pd.read_csv('/home/vikalin/Ira/+79.9657.csv', sep='\t')
df_1=ex_1.copy()
del df_1['localization score']
del df_1['spectrum']
df_1
mod_peptides=[]
peptides=[]
for elem in df_1['top isoform'].values:
    a=elem.find('[')
    if a!=-1:
        b=elem.find(']')
        peptide=elem[:a]+elem[b+1:]
    else:
        peptide=elem
    mod_peptide=str()   
    for i in range(len(peptide)):
        if peptide[i]=='T':
            mod_peptide=mod_peptide+peptide[i]+'*'
        else:
            mod_peptide=mod_peptide+peptide[i]
    mod_peptides.append(mod_peptide)
    peptides.append(peptide)
df_1['Peptide']=peptides        
df_1['Mod_Peptide']=mod_peptides
df_1['indexes']=peptides
count=[]
for elem in df_1['Mod_Peptide'].values:
    if '*' not in elem:
        count.append('0')
    else:
        count.append('1')
df_1['Count']=count        
df_1_1=df_1.where(df_1['Count']=='1')
df_1_1=df_1_1.dropna()
del df_1_1['top isoform']
del df_1_1['Count']
#Peptides=df_1_1.set_index('Peptide')
df_1_1['Protein']=None
Peptides=df_1_1.set_index('indexes')
#Peptides=df_1_1
Peptides


# In[ ]:


for description, sequence in fasta.read('/home/vikalin/Ira/sp_human_canonical_21032019_reversed.fasta'):
    print(description,sequence)


# In[115]:


#SRMS data
xl = pd.ExcelFile('/home/vikalin/Article/SRMS/Data.xlsx')
df1 = xl.parse('SRMS candidate substrates')
df1
#df1.columns=df1.loc[0]
#df=df1.drop(0)
#Peptides=((df.loc[df['Amb?'] == 'UNIQUE']).reset_index()).loc[:,df.columns.intersection(['Peptide']) ]
#indexes=[elem.replace('*','') for elem in Peptides['Peptide']]
#Peptides['Mod_Peptide']=indexes  
#Peptides.columns = ['Mod_Peptide','Peptide']
#Peptides['Protein']=None
#Peptides=Peptides.set_index('Peptide')
#Peptides


# In[3]:


df1['Leading proteins Uniprot ID']


# In[116]:


protein={}
for description, sequence in fasta.read(gzip.open('/home/vikalin/Article/HUMAN.fasta.gz','rt')):
    i=description.find('|')
    k=description.rfind('|')
    number=description[i+1:k]
    if number.find(';')==-1:   
        print(number)
        if number in df1['Leading proteins Uniprot ID'].values:
            protein[number]=sequence


# In[78]:


protein


# In[117]:


search=[]
for i in range(len(df1['Leading proteins Uniprot ID'].values)):
    if df1['Leading proteins Uniprot ID'].values[i] in protein:
        continue
    else:
        if (df1['Leading proteins Uniprot ID'].values[i]).count(';')==0:
            search.append(df1['Leading proteins Uniprot ID'].values[i])
print(search) 
print(len(search))


# In[118]:


import requests
def queryUniprot(accessionList, resFormat = "fasta"):
    """
    Retrieve sequences from uniprot by id list
    """
    params = {
        'from': 'ACC',
        'to': 'ACC',
        'format': resFormat,
        'query':' '.join(accessionList)
        }

    return requests.post('https://www.uniprot.org/uploadlists/', params)

a = queryUniprot(search)
print(a.text)


# In[119]:



b=a.text
search_1=b.split('>')
search_1=search_1[1:]
for elem in search_1:
    i=elem.find('|')+1
    k=elem.rfind('|')
    name=elem[i:k]
    l=elem.find('\n')
    seq=elem[l+1:]
    sequence=seq.replace('\n','')
    print(name,sequence)
    protein[name]=sequence


# In[120]:


proteins=[]
for i in range(len(df1['Leading proteins Uniprot ID'].values)):
    if df1['Leading proteins Uniprot ID'].values[i] in protein:
        proteins.append(protein[df1['Leading proteins Uniprot ID'].values[i]])
    else:
        proteins.append('-')
df1['Protein']=proteins
df1['Peptide']=proteins
df2=df1.where(df1['Protein']!='-')
df2=df2.dropna()
df2


# In[121]:


print(df2['Protein'].values[0])
a=df1['Protein'].values[0]
a[186]
mod_proteins=[]
for i in range(len(df2['Protein'].values)):
    a=df2['Protein'].values[i]
    b=df2['Position of phosphosite(s) on protein'].values[i]
    if type(b) == int:
        s=a[:b]+'*'+a[b:]
        print('s',s)
        mod_proteins.append(s)
    else:
        b=b.replace('.',',')
        b=b.replace(' ',',')
        B=b.split(',')
        B_1=[]
        for elem in B:
            if (elem.count('1')==0) and (elem.count('2')==0) and (elem.count('3')==0) and (elem.count('4')==0) and (
                elem.count('5')==0) and (elem.count('6')==0) and (elem.count('7')==0) and (elem.count('8')==0) and (
                elem.count('9')==0):
                continue
            else:
                B_1.append(int(elem))
        B_1.sort
        asterik=0
        s=a
        for elem in B_1:
            print('B',B_1,'elem',elem)
            number=int(elem)
            s=s[:number+asterik]+'*'+s[number+asterik:]
            asterik+=1
            print('s',s)
        mod_proteins.append(s)


# In[122]:


df2['Mod_Peptide']=mod_proteins
del df2['Protein names']
del df2['Gene names']
del df2['Position of phosphosite(s) on protein']
del df2['Amino acid carrying PTM']
del df2['Leading proteins Uniprot ID']
df2['index']=df2['Peptide'].values
df2=df2.set_index('index')
Peptides=df2
Peptides


# In[125]:


path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
#path_PTM='/home/vikalin/knowledgebase/PTM/PTM.txt'
path_PTM=None
print('Ищем модификации в:','\n','1-экспериментальном наборе','\n','2-базе ПТМ')
a=int(input())
if a==1:
    background_type='FASTA'
    print('Название исследуемого набора данных:')
    name_sample=str(input())
    print('Название организма:')
    organism=str(input())
    print('Введите длину интервала:')
    interval_length=int(input())
    print('Введите модификацию:')
    modification=str(input())
    print('Введите сайт модификации:')
    modification_site=str(input())
    chi2_results,chi2_selection,intervals,background=output_experimental(name_sample,Peptides,interval_length,modification_site,modification,path_FASTA,path_PTM)
elif a==2:
    background_type='FASTA'
    print('Название исследуемого набора данных:')
    name_sample=str(input())
    print('Название организма:')
    organism=str(input())
    print('Введите длину интервала:')
    interval_length=int(input())
    print('Введите модификацию:')
    modification=str(input())
    print('Введите сайт модификации:')
    modification_site=str(input())
    output_PTM(name_sample,organism,interval_length,modification_site,modification,path_FASTA,path_PTM)
    


# In[133]:


i=0
for elem in intervals:
    if elem[8]=='R':
        i+=1
print(i,len(intervals))        


# In[134]:


i=0
for elem in background:
    if elem[8]=='R':
        i+=1
print(i,len(background))


# In[139]:


file=open('/home/vikalin/Article/Ira_results_3/Phosphorylation/T_6_background:FASTA/intervals_T.txt','w')
for elem in intervals:
    file.write(elem+'\n')
file.close()    


# In[129]:


ex_1=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(хи-квадрат 18_06)/table6_quadruple.csv', sep=',')
double_proteom_FASTA=ex_1.copy()
double_motifs_FASTA=double_proteom_FASTA['Unnamed: 0'].values

ex_2=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(хи квадрат рандом 18_06)/table6_quadruple.csv', sep=',')
double_proteom_random=ex_2.copy()
double_motifs_random=double_proteom_random['Unnamed: 0'].values

ex_3=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(хи квадрат HeLa 18_06)/table6_quadruple.csv', sep=',')
double_proteom_HeLa=ex_3.copy()
double_motifs_HeLa=double_proteom_HeLa['Unnamed: 0'].values


# In[127]:


pip install matplotlib-venn


# In[130]:


#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
venn3([set(double_motifs_FASTA), set(double_motifs_random), set(double_motifs_HeLa)], set_labels = ('FASTA', 'random', 'HeLa'))
plt.title('Двухбуквенные мотивы \n')
#print('Повторяющиеся мотивы')
#for i in range(len(table_4)):
#    if (table_4[i] in table_2) and (table_4[i] in table_3):
#        print(table_4[i])table6_quadruple.csv
plt.savefig('/home/vikalin/double_chi.svg')


# In[ ]:


print(len(double_motifs_FASTA))
print(len(double_motifs_random))
print(len(double_motifs_HeLa))


# In[ ]:


#отцифровка таблицы с киназами
acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
n=[]
for i in range(len(acids)):
    a=[0]*(interval_length*2+1)
    n.append(a)
    
zero=['K','R']
for i in range(len(zero)):
    k=acids.index(zero[i])
    n[k][0]=1
    
one=['R','F','K','L','I','M','V','Y','E']
for i in range(len(one)):
    k=acids.index(one[i])
    n[k][1]=1

two=['R','F','K','L','N','H','V','Y','E','A']
for i in range(len(two)):
    k=acids.index(two[i])
    n[k][2]=1
    
three=['R','F','K','L','I','M','S','Y','E','T','A','Q','H','G','D']
for i in range(len(three)):
    k=acids.index(three[i])
    n[k][3]=1
    
four=['R','F','K','L','S','T','A','Q','E','H','G','D','P']
for i in range(len(four)):
    k=acids.index(four[i])
    n[k][4]=1
    
five=['R','G','K','L','I','M','N','S','E','T','A','F','D','Q']
for i in range(len(five)):
    k=acids.index(five[i])
    n[k][5]=1 
    
seven=['R','F','V','L','I','M','S','Y','E','P','Q','G']
for i in range(len(seven)):
    k=acids.index(seven[i])
    n[k][7]=1  
    
eight=['R','F','K','I','M','Y','E','V','G','D']
for i in range(len(eight)):
    k=acids.index(eight[i])
    n[k][8]=1    
    
nine=['R','F','K','L','I','M','S','E','Q','G','D','P','V']
for i in range(len(nine)):
    k=acids.index(nine[i])
    n[k][9]=1
    
ten=['R','F','K','L','I','M','S','Y','G','D','V','H']
for i in range(len(ten)):
    k=acids.index(ten[i])
    n[k][10]=1  
    
eleven=['F','L','A','P']
for i in range(len(eleven)):
    k=acids.index(eleven[i])
    n[k][11]=1    
    
path='/home/vikalin/'+'kinase'+'.txt'
saving=open(path,'w')
for i in range(0,21):
    for k in range(interval_length*2+1):
        saving.write(str(n[i][k])+' ')
    saving.write('\n')
saving.close()    


# In[ ]:


x=[]
for i in range(-interval_length,interval_length+1):
    if i==0:
        x.append(modification_site)
    else:    
        x.append(i)
fig = plt.figure()
ax = sns.heatmap(n,xticklabels=x,yticklabels=acids,linewidths=.5,cmap="RdBu_r",center=1,square=True)
fig.set_figwidth(15)    #  ширина и
fig.set_figheight(15)#  высота "Figure"
ax.set_title('Таблица киназ')

#сохранение рисунка

path='/home/vikalin/'+'kinase_new'+'.png'
plt.savefig(path)
plt.show()


# In[ ]:


comparison=[]
for i in range(len(acids)):
    a=[0]*(interval_length*2+1)
    comparison.append(a)
for i in range(0,21):
    for k in range(interval_length*2+1):
        if (n[i][k]==1) and (chi2_selection[i][k]==1):
            comparison[i][k]=50000
        elif (n[i][k]==1) and (chi2_selection[i][k]==0):
            comparison[i][k]=10000
        elif (n[i][k]==0) and (chi2_selection[i][k]==1):
            comparison[i][k]=5000    
            


# In[ ]:


x=[]
for i in range(-interval_length,interval_length+1):
    if i==0:
        x.append(modification_site)
    else:    
        x.append(i)
fig = plt.figure()
ax = sns.heatmap(comparison,xticklabels=x,yticklabels=acids,linewidths=.5,cmap="Paired",center=1,square=True)
fig.set_figwidth(15)    #  ширина и
fig.set_figheight(15)#  высота "Figure"
#ax.set_title(title)

#сохранение рисунка

path='/home/vikalin/'+'comparison_0.01'+'.png'
plt.savefig(path)
plt.show()


# In[92]:


kinases=dict()
kinases['CLK1']= [(13,0),(12,0),(13,2),(12,2),(13,4),(12,4),(12,9)] 
kinases['PKC-eta']=[(12,0),(12,3),(12,4),(13,4),(12,5),(2,7),(12,8),(12,9)]
kinases['PKC-alpha']=[(12,0),(12,1),(2,1),(12,3),(12,4),(13,4),(10,5),(2,7),(4,7),(12,8),(13,8),(12,9),(13,9),(16,9)]
kinases['PKA']=[(12,2),(13,2),(12,3),(13,3),(12,4),(13,4),(12,5),(17,5),(2,7),(5,7),(4,7),(6,7),(0,7),(5,8),(2,9),(19,10)]
kinases['Akt/PKB']=[(12,1),(12,3),(15,4),(14,4),(7,4),(15,5),(14,5),(7,5),(2,7),(4,7)]
kinases['PKC-delta']=[(12,1),(12,3),(13,4),(10,5),(2,7)]
kinases['PKC-gamma']=[(12,1),(12,2),(12,3),(13,4),(13,5),(10,5),(2,7),(12,8),(13,8),(12,9),(13,9),(13,10),(7,11)]
kinases['DMPK-E']=[(12,1),(13,1),(13,2),(12,3),(12,4),(12,5),(4,7),(6,7)]
kinases['Pim1']=[(12,1),(13,1),(12,2),(13,2),(12,3),(12,4),(13,4),(4,5)]
kinases['RSK1']=[(12,1),(13,1),(12,3),(12,4)]
kinases['PKC-epsilon']=[(13,2),(12,3),(12,4),(13,4),(16,4),(10,5),(6,7),(12,8),(12,9)]
kinases['SLK1']=[(12,2),(12,3),(2,4),(10,5),(2,7),(5,7),(4,7),(6,7),(0,7),(12,8),(12,9),(2,10),(5,10),(4,10),(6,10),(0,10)]
kinases['ZIPK']=[(12,2),(13,2),(12,3),(12,4),(12,5)]
kinases['NIMA']=[(17,2),(12,2),(2,3),(4,3),(3,3),(12,4),(13,4),(12,5),(13,5),(5,7),(6,7),(3,7),(12,7),(5,8),(6,8),(3,8),(12,8),(2,9),(5,9),(3,9),(6,9),(2,10),(5,10),(3,10)]
kinases['PKC-beta']=[(2,1),(4,1),(12,2),(13,2),(12,3),(13,4),(16,4),(10,5),(2,7),(3,7),(13,8),(13,9),(7,11)]
kinases['AMPK']=[(5,1),(4,1),(3,1),(6,1),(11,2),(13,2),(12,2),(11,3),(13,3),(12,3),(11,4),(13,4),(12,4),(5,10),(4,10),(3,10),(6,10)]
kinases['DCK1-b2']=[(2,1),(5,1),(4,1),(3,1),(6,1),(12,2),(12,3),(2,7),(5,7),(4,7),(3,7),(6,7)]
kinases['CHK1']=[(5,1),(4,1),(3,1),(6,1),(12,3),(13,3)]
kinases['CaMK1']=[(2,1),(5,1),(4,1),(3,1),(6,1),(12,3),(2,7),(5,7),(4,7),(3,7),(6,7)]
kinases['CaMK2']=[(2,1),(5,1),(4,1),(0,1),(6,1),(12,3),(13,3),(2,7),(5,7),(4,7),(3,7),(6,7),(19,8),(18,8)]
kinases['CaMK4']=[(2,1),(5,1),(4,1),(0,1),(6,1),(12,3)]
kinases['PKC-zeta']=[(2,1),(12,3),(2,7),(3,7),(2,8),(3,8)]
kinases['PKC-mu']=[(6,1),(4,1),(7,2),(4,2),(6,2),(12,3),(3,5)]
kinases['MSK1/2']=[(12,3)]
kinases['PAK']=[(12,3),(13,3),(12,4)]
kinases['Phos.Kinase']=[(12,3),(13,3),(2,5),(5,5),(3,5),(4,5),(2,7),(5,7),(3,7),(4,7),(6,7),(2,8),(13,8),(12,8),(5,9),(4,9),(2,10),(5,10),(4,10)]
kinases['PDK1']=[(2,2),(2,5),(2,7),(0,7)]
kinases['CK1d']=[(18,1),(18,2),(2,2),(19,3),(15,3),(7,4),(10,4),(14,4),(10,5),(5,7),(2,8),(10,8),(5,8),(0,8),(10,9),(2,9),(5,9),(10,10),(2,10),(2,11),(4,11),(9,11)]
kinases['CK1g']=[(0,1),(18,2),(0,2),(0,3),(19,3),(15,3),(7,4),(19,4),(7,5),(10,5),(5,7),(2,8),(10,8),(5,8),(0,8),(10,9),(2,9),(5,9),(10,10),(2,10),(2,11),(4,11),(9,11)]
kinases['BARK']=[(18,4)]
kinases['ATM']=[(5,4),(4,4),(3,4),(9,4),(18,5),(19,5),(16,7),(18,8)]
kinases['DNAPK']=[(18,5),(16,7),(18,8)]
kinases['CK2']=[(19,8),(18,8),(19,9),(18,9),(15,9),(19,10),(18,10)]
kinases['ERK1']=[(9,4),(9,7),(2,10),(5,10),(0,10)]
kinases['p38 MAPK']=[(9,4),(9,7),(10,3),(16,5),(3,5),(5,8)]
kinases['CDK2']=[(9,4),(9,7),(12,2),(3,5),(12,9),(13,9),(13,10)]
kinases['CDK4']=[(9,4),(9,7),(4,5),(9,9),(13,10),(12,10),(11,10)]
kinases['CDK5']=[(9,4),(9,7),(13,9)]
kinases['CDK1']=[(12,5),(13,5),(9,7),(12,8),(13,8),(12,9),(13,9)]
kinases['GSK3']=[(15,10),(9,11)]
kinases['TGF-betaR1']=[(15,7),(15,9)]

#только для pT
# kinases['MAPKAPK2']=[(5,1),(4,1),(6,1),(12,3),(4,5),(5,7),(4,7),(3,7),(6,7)]
# kinases['LKB1']=[(4,4)]
# kinases['mTOR']=[(2,5),(0,7)]
# kinases['MEK3']=[(10,7),(0,8)]
# kinases['SIM']=[(18,7),(0,8)]


# In[93]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы+хи квадрат)/table6_double.csv', sep=',')
double_proteom=ex.copy()
a=list(kinases.keys())
table_2_chi2=dict()
print(len(double_proteom['Location'].values))
#фиксируем конкретный мотив
for i in range(len(double_proteom['Location'].values)):
    s=(double_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    acid_1=(acid_1_1,acid_1_2)
    #print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-2])
    acid_2=(acid_2_1,acid_2_2)
    #бежим по киназам
    for elem in a:
        if (acid_1 in kinases[elem]) and (acid_2 in kinases[elem]):
            if elem in table_2_chi2:
                table_2_chi2[elem]=table_2_chi2[elem]+' '+double_proteom['Motifs'].values[i]
            else:
                table_2_chi2[elem]=double_proteom['Motifs'].values[i]
table_2_chi2                
#del double_proteom['index']


# In[94]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы+хи квадрат)/table6_triple.csv', sep=',')
double_proteom=ex.copy()
a=list(kinases.keys())
table_3_chi2=dict()
print(len(double_proteom['Location'].values))
#фиксируем конкретный мотив
for i in range(len(double_proteom['Location'].values)):
    s=(double_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    acid_1=(acid_1_1,acid_1_2)
    #print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-1])
    acid_2=(acid_2_1,acid_2_2)
    #print(acid_2)
    acid_3_1=int(s[4][2:])
    acid_3_2=int(s[5][:-2])
    acid_3=(acid_3_1,acid_3_2)
    #print(acid_3)
    #бежим по киназам
    for elem in a:
        if (acid_1 in kinases[elem]) and (acid_2 in kinases[elem]) and (acid_3 in kinases[elem]):
            if elem in table_3_chi2:
                table_3_chi2[elem]=table_3_chi2[elem]+' '+double_proteom['Motifs'].values[i]
            else:
                table_3_chi2[elem]=double_proteom['Motifs'].values[i]
table_3_chi2 


# In[95]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы+хи квадрат)/table6_quadruple.csv', sep=',')
double_proteom=ex.copy()
a=list(kinases.keys())
table_4_chi2=dict()
print(len(double_proteom['Location'].values))
#фиксируем конкретный мотив
for i in range(len(double_proteom['Location'].values)):
    s=(double_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    acid_1=(acid_1_1,acid_1_2)
    #print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-1])
    acid_2=(acid_2_1,acid_2_2)
    #print(acid_2)
    acid_3_1=int(s[4][2:])
    acid_3_2=int(s[5][:-1])
    acid_3=(acid_3_1,acid_3_2)
    #print(acid_3)
    acid_4_1=int(s[6][2:])
    acid_4_2=int(s[7][:-2])
    acid_4=(acid_4_1,acid_4_2)
    #print(acid_4)
    #бежим по киназам

    for elem in a:
        if (acid_1 in kinases[elem]) and (acid_2 in kinases[elem]) and (acid_3 in kinases[elem]) and (acid_4 in kinases[elem]):
            print('!')
            if elem in table_4_chi2:
                table_4_chi2[elem]=table_4_chi2[elem]+' '+double_proteom['Motifs'].values[i]
            else:
                table_4_chi2[elem]=double_proteom['Motifs'].values[i]
table_4_chi2 


# In[96]:


table_chi2=dict()
for elem in list(table_2_chi2.keys()):
    table_chi2[elem]=table_2_chi2[elem]
for elem in list(table_3_chi2.keys()):
    if elem not in table_chi2:
        table_chi2[elem]=table_3_chi2[elem]
    else:
        table_chi2[elem]=table_chi2[elem]+' '+table_3_chi2[elem]
table_chi2    


# In[97]:



motifs=[]
for elem in list(kinases.keys()):
    if elem in list(table_chi2.keys()):
        motif=(table_chi2[elem]).replace(' ',',')
        motifs.append(motif)
    else:
        motifs.append(None)
TABLE=pd.DataFrame({'Kinase':list(kinases.keys()),'Motifs_chi2':motifs})
TABLE


# In[98]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы + bi)/table6_double.csv', sep=',')
double_proteom=ex.copy()
double_proteom
count=[]
motifs=[]
for elem in double_proteom['Unnamed: 0'].values:
    if elem not in motifs:
        count.append(0)
        motifs.append(elem)
    else:
        count.append(1)
double_proteom['Count']=count
double_proteom=double_proteom[double_proteom['Count']==0]
del double_proteom['Count']
#double_proteom=double_proteom.where(double_proteom['Probability']<0.05/493)
#double_proteom=double_proteom.dropna()
double_proteom


# In[99]:


a=list(kinases.keys())
table_2_bi=dict()
#фиксируем конкретный мотив
for i in range(len(double_proteom['Location'].values)):
    s=(double_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    acid_1=(acid_1_1,acid_1_2)
    #print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-2])
    acid_2=(acid_2_1,acid_2_2)
    #бежим по киназам
    for elem in a:
        if (acid_1 in kinases[elem]) and (acid_2 in kinases[elem]):
            if elem in table_2_bi:
                table_2_bi[elem]=table_2_bi[elem]+' '+double_proteom['Unnamed: 0'].values[i]
            else:
                table_2_bi[elem]=double_proteom['Unnamed: 0'].values[i]
table_2_bi


# In[100]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы + bi)/table6_triple.csv', sep=',')
triple_proteom=ex.copy()
triple_proteom
count=[]
motifs=[]
for elem in triple_proteom['Unnamed: 0'].values:
    if elem not in motifs:
        count.append(0)
        motifs.append(elem)
    else:
        count.append(1)
triple_proteom['Count']=count
triple_proteom=triple_proteom[triple_proteom['Count']==0]
del triple_proteom['Count']
#triple_proteom=triple_proteom.where(triple_proteom['Occurrences']>=10)
#triple_proteom=triple_proteom.dropna()
#triple_proteom=triple_proteom.where(triple_proteom['Probability']<0.05/308)
#triple_proteom=triple_proteom.dropna()
triple_proteom


# In[101]:


table_3_bi=dict()
#фиксируем конкретный мотив
for i in range(len(triple_proteom['Location'].values)):
    comb=[]
    s=(triple_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    if acid_1_2!=6:    
        acid_1=(acid_1_1,acid_1_2)
        comb.append(acid_1)
    #print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-1])
    if acid_2_2!=6:   
        acid_2=(acid_2_1,acid_2_2)
        comb.append(acid_2)
    #print(acid_2)
    acid_3_1=int(s[4][2:])
    acid_3_2=int(s[5][:-1])
    if acid_3_2!=6:   
        acid_3=(acid_3_1,acid_3_2)
        comb.append(acid_3)
    #print(acid_3)
    acid_4_1=int(s[6][2:])
    acid_4_2=int(s[7][:-2])
    if acid_4_2!=6:   
        acid_4=(acid_4_1,acid_4_2)
        comb.append(acid_4)
    #print(acid_3)
    
    #бежим по киназам
    for elem in a:
        if (comb[0] in kinases[elem]) and (comb[1] in kinases[elem]) and (comb[2] in kinases[elem]):
            print(comb)
            if elem in table_3_bi:
                table_3_bi[elem]=table_3_bi[elem]+' '+triple_proteom['Unnamed: 0'].values[i]
            else:
                table_3_bi[elem]=triple_proteom['Unnamed: 0'].values[i]
table_3_bi 


# In[102]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы + bi)/table6_quadruple.csv', sep=',')
triple_proteom=ex.copy()
triple_proteom
count=[]
motifs=[]
for elem in triple_proteom['Unnamed: 0'].values:
    if elem not in motifs:
        count.append(0)
        motifs.append(elem)
    else:
        count.append(1)
triple_proteom['Count']=count
triple_proteom=triple_proteom[triple_proteom['Count']==0]
del triple_proteom['Count']
triple_proteom
table_4_bi=dict()
#фиксируем конкретный мотив
for i in range(len(triple_proteom['Location'].values)):
    comb=[]
    s=(triple_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    if acid_1_2!=6:    
        acid_1=(acid_1_1,acid_1_2)
        comb.append(acid_1)
    print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-1])
    if acid_2_2!=6:   
        acid_2=(acid_2_1,acid_2_2)
        comb.append(acid_2)
    print(acid_2)
    acid_3_1=int(s[4][2:])
    acid_3_2=int(s[5][:-1])
    if acid_3_2!=6:   
        acid_3=(acid_3_1,acid_3_2)
        comb.append(acid_3)
    print(acid_3)
    acid_4_1=int(s[6][2:])
    acid_4_2=int(s[7][:-1])
    if acid_4_2!=6:   
        acid_4=(acid_4_1,acid_4_2)
        comb.append(acid_4)
    print(acid_4)
    acid_5_1=int(s[8][2:])
    acid_5_2=int(s[9][:-2])
    if acid_5_2!=6:   
        acid_5=(acid_5_1,acid_5_2)
        comb.append(acid_5)
    print(acid_5)
    
    #бежим по киназам
    for elem in a:
        if (comb[0] in kinases[elem]) and (comb[1] in kinases[elem]) and (comb[2] in kinases[elem]) and (comb[3] in kinases[elem]):
            print('!')
            if elem in table_4_bi:
                table_4_bi[elem]=table_4_bi[elem]+' '+triple_proteom['Unnamed: 0'].values[i]
            else:
                table_4_bi[elem]=triple_proteom['Unnamed: 0'].values[i]
table_4_bi 


# In[103]:


table_bi=dict()
for elem in list(table_2_bi.keys()):
    table_bi[elem]=table_2_bi[elem]
for elem in list(table_3_bi.keys()):
    if elem not in table_bi:
        table_bi[elem]=table_3_bi[elem]
    else:
        table_bi[elem]=table_bi[elem]+' '+table_3_bi[elem]
table_bi


# In[104]:


motifs=[]
for elem in list(kinases.keys()):
    if elem in list(table_bi.keys()):
        motif=(table_bi[elem]).replace(' ',',')
        motifs.append(motif)
    else:
        motifs.append(None)


# In[105]:



TABLE['Motifs_bi']=motifs
TABLE.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/comparison_chi2_bi_T.csv')
TABLE.to_excel('/home/vikalin/Article/example/PHOSPHORYLATION/comparison_chi2_bi_T.xlsx')
TABLE


# In[107]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
names=[]
for i in range(len(TABLE['Kinase'].values)):
    a=TABLE['Kinase'].values[i]
    b=kinases[a]
#    print(b)
    kinasa=dict()
    for elem in b:
#        print(elem)
        acid,pos=elem
        if pos in kinasa:
#            print(kinasa[pos])
            kinasa[pos]=kinasa[pos]+','+acids[acid]
        else:
            kinasa[pos]=acids[acid]
    print(kinasa)        
     #собираем мотив
    seq=''
    for i in range(13):
        if i in kinasa:
            motifs=kinasa[i]
            motifs=motifs.replace(',','')
            print('motifs',motifs)
            if len(motifs)>1:
                seq=seq+'['+motifs+']'
            else:
                seq=seq+motifs
        else:
            if i==6:
                seq=seq+'s'
            else:
                seq=seq+'.'
    print(seq)    
            
    names.append(seq)            


# In[108]:


TABLE_1=TABLE.copy()
count=[]
TABLE_1['Motif']=names
for i in range(len(TABLE_1['Kinase'].values)):
    if TABLE_1['Motifs_chi2'].values[i]==None and TABLE_1['Motifs_bi'].values[i]==None:
        count.append(0)
    else:
        count.append(1)
TABLE_1['Count']=count
for i in range(len(TABLE_1['Kinase'].values)):
    if TABLE_1['Motifs_chi2'].values[i]==None:
        TABLE_1['Motifs_chi2'].values[i]='-'
    if TABLE_1['Motifs_bi'].values[i]==None:
        TABLE_1['Motifs_bi'].values[i]='-'   
        
TABLE_2=TABLE_1.where(TABLE_1['Count']==1).dropna()
del TABLE_2['Count']
TABLE_2 = TABLE_2[['Kinase','Motif', 'Motifs_chi2','Motifs_bi']]
TABLE_2=TABLE_2.reset_index()
del TABLE_2['index']
TABLE_2.rename(columns={'Kinase':'Киназы','Motif':'Полный мотив','Motifs_chi2':'Мотив хи-квадрат','Motifs_bi':'Мотивы бином'}, inplace=True)
TABLE_2


# In[21]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/comparison_chi2_bi_T.csv', sep=',')
double_proteom=ex.copy()
del double_proteom['Unnamed: 0']
count=[]
motifs_bi=[]
for i in range(len(double_proteom['Motifs_chi2'].values)):
    if (type(double_proteom['Motifs_chi2'].values[i]) is float) and (type(double_proteom['Motifs_bi'].values[i]) is float):
        count.append(1)
        motifs_bi.append('-')
    else:
        count.append(0)
        if (type(double_proteom['Motifs_bi'].values[i]) is float):
            motifs_bi.append('-')
        else:
            motifs_bi.append(double_proteom['Motifs_bi'].values[i])
double_proteom['Count']=count
del double_proteom['Motifs_bi']
double_proteom['Motifs_bi']=motifs_bi
double_proteom_1=double_proteom.where(double_proteom['Count']==0).dropna()
del double_proteom_1['Count']
double_proteom_1
#motifs_bi
double_proteom_1=double_proteom_1.reset_index()
del double_proteom_1['index']
double_proteom_1.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/comparison_chi2_bi_final_T.csv')
double_proteom_1.rename(columns={'Kinase':'Киназы','Motifs_chi2':'Мотивы хи-квадрат','Motifs_bi':'Мотивы бином'})


# In[41]:


print(len(motifs_bi))
#TABLE['Motifs_bi']=motifs_bi
#TABLE


# In[7]:


#проверка второго экспериментального датасета
exp_2=['K.y','K...y','yG','yA','R.y','R...y']
ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/Y_10_background:FASTA(хи квадрат)/table10_primary.csv', sep=',')
chi_proteom=ex.copy()
chi=[]
for i in range(len(chi_proteom['Motifs'].values)):
    if chi_proteom['Motifs'].values[i] in exp_2:
        chi.append(chi_proteom['Motifs'].values[i])
ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/Y_10_background:FASTA(bi)/table10_single.csv', sep=',')
bi_proteom=ex.copy()
bi_proteom
bi=[]
for i in range(len(bi_proteom['Unnamed: 0'].values)):
    if bi_proteom['Unnamed: 0'].values[i] in exp_2:
        bi.append(bi_proteom['Unnamed: 0'].values[i])
table=pd.DataFrame({'Киназа':'PTK70','Мотивы хи-квадрат':chi,'Мотивы бином':bi})  
table


# In[46]:


double_proteom['Motifs'].values


# In[47]:


a=list(kinases.keys())
table_2_chi2=dict()
#фиксируем конкретный мотив
for i in range(len(double_proteom['Location'].values)):
    s=(double_proteom['Location'].values[i]).split(',')
    acid_1_1=int(s[0][2:])
    acid_1_2=int(s[1][:-1])
    acid_1=(acid_1_1,acid_1_2)
    #print(acid_1)
    acid_2_1=int(s[2][2:])
    acid_2_2=int(s[3][:-2])
    acid_2=(acid_2_1,acid_2_2)
    #бежим по киназам
    for elem in a:
        if (acid_1 in kinases[elem]) and (acid_2 in kinases[elem]):
            print(double_proteom['Motifs'].values[i])
            if elem in table_2_chi2:
                table_2_chi2[elem]=table_2_chi2[elem]+' '+double_proteom['Motifs'].values[i]
            else:
                table_2_chi2[elem]=double_proteom['Motifs'].values[i]
                


# In[44]:


table_2_chi2


# In[12]:


len(double_proteom['Location'].values[i])


# In[289]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
P_HeLa=[]
acids_HeLa=[]
occurrences_HeLa=[]
location_HeLa=[]
for i in range(len(acids)):
    for k in range(13):
#        acids_HeLa.append(acids[i]+str((i,k)))
#        P_HeLa.append(P[i][k])
#        occurrences_HeLa.append(occurrences[i][k])
        location_HeLa.append((i,k))


# In[540]:


primary_HeLa=pd.DataFrame({'Acid':acids_HeLa,'P':P_HeLa,'Occurrences':occurrences_HeLa})
sorted_primary_HeLa=primary_HeLa.copy()
sorted_primary_HeLa.sort_values('P',ascending=True,inplace=True)
sorted_primary_HeLa


# In[8]:


file=open('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA (P, occurrences)/P.txt','r')
P_FASTA=[]
for line in file:
    t=line.split(' ')
    for i in range(len(t)):
        if t[i]!='\n':
            P_FASTA.append(float(t[i]))        


# In[55]:


#РАБОТАЮ СЕЙЧАС ЗДЕСЬ
file=open('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA/occurrences.txt','r')
occurrences_FASTA=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    occurrences_FASTA.append(a)
k=0    
for line in file:
    t=line.split(' ')
    for i in range(len(t)):
        if t[i]!='\n':
            occurrences_FASTA[k][i]=int(t[i])
    k+=1        
occurrences_FASTA


# In[56]:


all_exp=0
for k in range(0,21):
    all_exp+=occurrences_FASTA[k][0]
all_exp    


# In[57]:


file=open('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA/background.txt','r')
background_FASTA=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    background_FASTA.append(a)
k=0    
for line in file:
    t=line.split(' ')
    for i in range(len(t)):
        if t[i]!='\n':
            background_FASTA[k][i]=int(t[i])
    k+=1        
background_FASTA


# In[58]:


all_back=0
for k in range(0,21):
    all_back+=background_FASTA[k][0]
all_back


# In[59]:


norm_FASTA=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    norm_FASTA.append(a)
for i in range(0,21):
    for k in range(0,13):
        norm_FASTA[i][k]=(background_FASTA[i][k])/all_back


# In[60]:


expected_FASTA=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    expected_FASTA.append(a)
for i in range(0,21):
    for k in range(0,13):
        expected_FASTA[i][k]=(norm_FASTA[i][k])*all_exp


# In[67]:


from scipy.stats import chisquare
chisquare([occurrences_FASTA[9][0],all_exp-occurrences_FASTA[9][0]],f_exp=[expected_FASTA[9][0],all_exp-expected_FASTA[9][0]])


# In[144]:


chi2_results=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    chi2_results.append(a)
for i in range(0,21):
    for k in range(0,13):
        if k!=6:
            chisq,p_value=chisquare([occurrences_FASTA[i][k],all_exp-occurrences_FASTA[i][k]],f_exp=[expected_FASTA[i][k],all_exp-expected_FASTA[i][k]])
            chi2_results[i][k]=p_value


# In[155]:


chi2_selection=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    chi2_selection.append(a)
for i in range(0,21):
    for k in range(0,13):
        if (chi2_results[i][k]<0.01/(21*13)) and (k!=6) and (occurrences_FASTA[i][k]>10):
            chi2_selection[i][k]=1


# In[156]:


chi2_selection


# In[157]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
x=[]
for i in range(-6,6+1):
    if i==0:
        x.append('S')
    else:    
        x.append(i)
fig = plt.figure()
ax = sns.heatmap(chi2_selection,xticklabels=x,yticklabels=acids,linewidths=.5,cmap="RdBu_r",center=1,square=True)
fig.set_figwidth(15)    #  ширина и
fig.set_figheight(15)#  высота "Figure"

path='/home/vikalin/result_chi2_0.01'+'.png'
plt.savefig(path)
plt.show()


# In[302]:


expected_FASTA


# In[304]:


hi_q=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    hi_q.append(a)
for i in range(0,21):
    for k in range(0,13):
        if k!=6 and expected_FASTA[i][k]!=0:
            hi_q[i][k]=(((occurrences_FASTA[i][k]-expected_FASTA[i][k])**2)/expected_FASTA[i][k])+(((all_exp-occurrences_FASTA[i][k]-(all_exp-expected_FASTA[i][k]))**2)/(all_exp-expected_FASTA[i][k]))


# In[305]:


hi_q


# In[353]:


hi_q_threshold=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    hi_q_threshold.append(a)
for i in range(0,21):
    for k in range(0,13):
        if (k!=6) and (hi_q[i][k]>40): #and (occurrences_FASTA[i][k]>10):
            hi_q_threshold[i][k]=1
            
            


# In[354]:


hi_q_threshold


# In[355]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
x=[]
for i in range(-6,6+1):
    if i==0:
        x.append('S')
    else:    
        x.append(i)
fig = plt.figure()
ax = sns.heatmap(hi_q_threshold,xticklabels=x,yticklabels=acids,linewidths=.5,cmap="RdBu_r",center=1,square=True)
fig.set_figwidth(15)    #  ширина и
fig.set_figheight(15)#  высота "Figure"
#ax.set_title(title)

#сохранение рисунка
#if name!=None:
#    path=path_results+'/'+name+'_'+modification_site+str(interval_length)+'.png'
#    plt.savefig(path)
plt.show()


# In[340]:


hi_q


# In[3]:


print(occurrences_FASTA)


# In[ ]:





# In[4]:


print(len(occurrences_FASTA))


# In[5]:


type(P_FASTA[0])


# In[12]:


#primary_FASTA=pd.DataFrame({'Acid':acids_HeLa,'P':P_FASTA,'Occurrences':occurrences_FASTA})
primary_FASTA=pd.DataFrame({'Location':location_HeLa,'P':P_FASTA,'Occurrences':occurrences_FASTA})
sorted_primary_FASTA=primary_FASTA.copy()
sorted_primary_FASTA.sort_values('P',ascending=True,inplace=True)
sorted_primary_FASTA


# In[33]:


count=[]
for elem in sorted_primary_FASTA['Location'].values:
    i,k=elem
    if (i==20) or (k==6):
        count.append(1)
    else:
        count.append(0)
sorted_primary_FASTA['Count']=count
sorted_primary_FASTA_1=sorted_primary_FASTA[sorted_primary_FASTA['Count']==0]
sorted_primary_FASTA_1


# In[34]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
interval_length=6
modification_site='S'
indexes=[]
for elem in sorted_primary_FASTA_1['Location'].values:
    i,j=elem
    if j<interval_length:
        s = [acids[i], '.'*(interval_length-j-1), modification_site.lower()]
        index=''.join(s)
    #index=acids[i]+'.'*(interval_length-j-1)+modification_site.lower()
    else:
        s = [modification_site.lower(), '.'*(j-interval_length-1), acids[i]]
        index=''.join(s)
    indexes.append(index)    


# In[35]:


sorted_primary_FASTA_1['Motifs']=indexes
sorted_primary_FASTA_1


# In[44]:


scores=[]
for elem in sorted_primary_FASTA_1['P'].values:
    if elem==0:
        score='P=0'
    else:    
        score=-math.log10(elem)
    scores.append(score)
sorted_primary_FASTA_1['Score_FASTA']=scores
del sorted_primary_FASTA_1['Count']


# In[46]:


del sorted_primary_FASTA_1['P']


# In[14]:


print(len(indexes))


# In[380]:


print(sorted_primary_FASTA)


# In[544]:


Probabilities=[]
Matches=[]
q=-40
while q<=0:
    P=10**(q)
    sorted_primary_HeLa_1=sorted_primary_HeLa.where(sorted_primary_HeLa['P']<=P)
    sorted_primary_HeLa_1=sorted_primary_HeLa_1.dropna()
    print(sorted_primary_HeLa_1)
    sorted_primary_FASTA_1=sorted_primary_FASTA.where(sorted_primary_FASTA['P']<=P)
    sorted_primary_FASTA_1=sorted_primary_FASTA_1.dropna()
    print(sorted_primary_FASTA_1)
    i=0
    if len(sorted_primary_FASTA_1['Acid'].values)!=0:
        for elem in sorted_primary_FASTA_1['Acid'].values:
            if (elem not in sorted_primary_HeLa_1['Acid'].values):
                i+=1
                print(elem)
    Probabilities.append(math.log10(P))
    Matches.append(i)
    q=q+1 
graph=pd.DataFrame({'Differences':Matches}, index=Probabilities)
graph.plot()
plt.show() 


# In[483]:


Probabilities


# In[77]:


#считаем количество интервалов с данным комплексным мотивом в background
x=sum(1 for interval in background if ((interval[7]=='P') and (interval[9]=='R')))
print('x:',x)

#считаем p-value этого комплексного мотива
p=x/len(background)
print('p:',p)

#теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
n=sum(1 for interval in intervals if ((interval[7]=='P') and (interval[9]=='R')))
print('n:',n)

#считаем Р вероятность по биномиальной формуле для комплексного мотива
P_AB=0
c=n
while c<=len(intervals):
    P_AB=P_AB+binom.pmf(n,len(intervals),p,loc=0)
    c+=1
print('P_AB:',P_AB)    
#print('P_AB:',P_AB)

#нашли вероятность P(AB), теперь найдем условную вероятность
P_B=P[12][9]
print('P_B:',P_B)
#print('P_B:',P_B,motif_1_x,motif_1_y)
P_A_B=P_AB/P_B
print('P_A_B:',P_A_B)


# In[126]:


#считаем количество интервалов с тройным мотивом в background
x=sum(1 for interval in background if ((interval[5]=='D') 
                                       and (interval[7]=='E') 
                                       and (interval[9]=='E')))
print('x:',x)
#считаем p-value этого комплексного мотива
p=x/len(background)
print('p:',p)

#теперь нужно подсчитать встречаемость комплексного мотива в исходном датасете
n=sum(1 for interval in intervals if ((interval[5]=='D') 
                                       and (interval[7]=='E') 
                                       and (interval[9]=='E')))
print('n:',n)

#считаем Р вероятность по биномиальной формуле для тройного мотива
P_ABC=0
c=n
while c<=len(intervals):
    P_ABC=P_ABC+binom.pmf(n,len(intervals),p,loc=0)
    c+=1
print('P_ABC:',P_ABC)    
#print('P_AB:',P_AB)

#нашли вероятность P(ABC), теперь найдем условную вероятность
P_BC=P_AB
print('P_BC:',P_BC)
#print('P_B:',P_B,motif_1_x,motif_1_y)
P_A_BC=P_ABC/P_BC
print('P_A_BC:',P_A_BC)
#print('P_A_B:',P_A_B)


# In[136]:


pip install matplotlib-venn


# In[137]:


pip install matplotlib-venn
#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[138]:


xl = pd.ExcelFile('/home/vikalin/Article/example/table 1.xls')
df1 = xl.parse('Sheet1')
df1.columns=df1.loc[2]
df1=df1[3:24]
df1
df=df1['Motif']
df_n=df.T.groupby(level=0).first().T
table_1=df_n['Motif'].values
table_1


# In[139]:


xl = pd.ExcelFile('/home/vikalin/Article/example/table 2.xls')
df1 = xl.parse('Sheet1')
df1.columns=df1.loc[2]
df1
df1=df1[3:21]
df1
df=df1['Motif']
df_n=df.T.groupby(level=0).first().T
table_2=df_n['Motif'].values
table_2


# In[140]:


xl = pd.ExcelFile('/home/vikalin/Article/example/table 3.xls')
df1 = xl.parse('Sheet1')
df1.columns=df1.loc[2]
df1
df1=df1[3:20]
df1
df=df1['Motif']
df_n=df.T.groupby(level=0).first().T
table_3=df_n['Motif'].values
table_3


# In[141]:


xl = pd.ExcelFile('/home/vikalin/Article/example/table 4.xls')
df1 = xl.parse('Sheet1')
df1.columns=df1.loc[2]
df1
df1=df1[3:20]
df1
df=df1['Motif']
df_n=df.T.groupby(level=0).first().T
table_4=df_n['Motif'].values
table_4


# In[148]:


table=[]


# In[152]:


for i in range(len(table_4)):
    if table_4[i] not in table:
        table.append(table_4[i])


# In[ ]:





# In[ ]:





# In[91]:


pip install matplotlib-venn
#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
venn3([set(table_4), set(table_2), set(table_3)], set_labels = ('Half sized background database', 'Proteomic background database', 'Random background database'))
plt.title('Мотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(table_4)):
    if (table_4[i] in table_2) and (table_4[i] in table_3):
        print(table_4[i])


# In[ ]:





# In[84]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(random)/table6_single.csv', sep=',')
single_random=ex.copy()
single_random.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
#motifs_single_Hela=single_HeLa['Motifs'].values
single_random


# In[3]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(P,occurrences,background)/table6_single.csv', sep=',')
single_proteom=ex.copy()
single_proteom.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
single_proteom
#motifs_single_proteom=single_proteom['Motifs'].values


# In[4]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(P,occurrences,background)/table6_double.csv', sep=',')
double_proteom=ex.copy()
double_proteom.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
double_proteom


# In[ ]:


for i in range(len(double)


# In[10]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(random)/table6_single.csv', sep=',')
single_random=ex.copy()
single_random.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
motifs_single_random=single_random['Motifs'].values


# In[13]:


venn3([set(motifs_single_Hela), set(motifs_single_proteom), set(motifs_single_random)], set_labels = ('Background:HeLa', 'Background:Proteom', 'Background:Random'))
plt.title('Первичные мотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_single_Hela)):
    if (motifs_single_Hela[i] in motifs_single_proteom) and (motifs_single_Hela[i] in motifs_single_random):
        print(motifs_single_Hela[i])


# In[14]:


venn2([set(motifs_single_Hela), set(motifs_single_proteom)], set_labels = ('Background:HeLa', 'Background:Proteom'))
plt.title('Первичные мотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_single_Hela)):
    if (motifs_single_Hela[i] in motifs_single_proteom):
        print(motifs_single_Hela[i])


# In[15]:


venn2([set(motifs_single_Hela), set(motifs_single_random)], set_labels = ('Background:HeLa', 'Background:Random'))
plt.title('Первичные мотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_single_Hela)):
    if (motifs_single_Hela[i] in motifs_single_random):
        print(motifs_single_Hela[i])


# In[16]:


venn2([set(motifs_single_proteom), set(motifs_single_random)], set_labels = ( 'Background:Proteom', 'Background:Random'))
plt.title('Первичные мотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_single_proteom)):
    if (motifs_single_proteom[i] in motifs_single_random):
        print(motifs_single_proteom[i])


# In[88]:



ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:HeLa/table6_triple.csv', sep=',')
triple_HeLa=ex.copy()
triple_HeLa.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
del triple_HeLa['P_ABC']
del triple_HeLa['P_BC']
#double_HeLa=double_HeLa[:150]
#print(len(set(double_HeLa['Motifs'].values)))
new_column=[]
for elem in triple_HeLa['Probability'].values:
    if elem!=math.inf:       
        i=round(-math.log(elem))
    else:
        i=0
    new_column.append(i)
triple_HeLa['Score']=new_column
del triple_HeLa['Probability']
cols = list(triple_HeLa.columns)
a, b = cols.index('Score'), cols.index('Location')
cols[b], cols[a] = cols[a], cols[b]
triple_HeLa = triple_HeLa[cols]
i=[]
counter=[]
a=triple_HeLa['Motifs'].values
for elem in a:
    if elem not in i:
        i.append(elem)
        counter.append(1)
    else:
        counter.append(0)
triple_HeLa['Counter']=counter
triple_HeLa=triple_HeLa[triple_HeLa['Counter']==1]
del triple_HeLa['Counter']   
#double_HeLa=double_HeLa[:40]
triple_HeLa
triple_HeLa.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/triple_HeLa_all.csv')
#motifs_triple_HeLa=triple_HeLa['Motifs'].values


# In[90]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA (proteom)/table6_triple.csv', sep=',')
triple_proteom=ex.copy()
triple_proteom.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
del triple_proteom['P_ABC']
del triple_proteom['P_BC']
#double_proteom=double_proteom[:150]
#print(len(set(double_proteom['Motifs'].values)))
new_column=[]
for elem in triple_proteom['Probability'].values:
    if elem!=math.inf:  
        i=round(-math.log(elem))
    else:
        i=0
    new_column.append(i)
triple_proteom['Score']=new_column
del triple_proteom['Probability']
cols = list(triple_proteom.columns)
a, b = cols.index('Score'), cols.index('Location')
cols[b], cols[a] = cols[a], cols[b]
triple_proteom = triple_proteom[cols]
i=[]
counter=[]
a=triple_proteom['Motifs'].values
for elem in a:
    if elem not in i:
        i.append(elem)
        counter.append(1)
    else:
        counter.append(0)
triple_proteom['Counter']=counter
triple_proteom=triple_proteom[triple_proteom['Counter']==1]
del triple_proteom['Counter']   
#triple_proteom=triple_proteom[:40]
triple_proteom.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/triple_proteom_all.csv')
triple_proteom
#motifs_triple_proteom=triple_proteom['Motifs'].values


# In[73]:


Probabilities=[]
Matches=[]
q=130
P=10**(-q)
while q>=-20:
    P=10**(-q)
    double_HeLa_exp=double_HeLa.where(double_HeLa['Probability']<=P)
    double_HeLa_exp=double_HeLa_exp.dropna()
    double_proteom_exp=double_proteom.where(double_proteom['Probability']<=P)
    double_proteom_exp=double_proteom_exp.dropna()
    i=0
    if len(double_proteom_exp['Motifs'].values)!=0:
        for elem in double_proteom_exp['Motifs'].values:
            if (elem!=None) and (elem not in double_HeLa_exp['Motifs'].values):
                i+=1
    Probabilities.append(math.log(P))
    Matches.append(i)
    q=q-10 
graph=pd.DataFrame({'Differences':Matches}, index=Probabilities)
graph.plot()
plt.show()


# In[74]:


Matches


# In[69]:


double_proteom_exp=double_proteom.where(double_proteom['Probability']<=10**(-130))
#double_proteom_exp=double_proteom.where(double_proteom['Probability']<=P)
double_proteom_exp=double_proteom_exp.dropna()
double_proteom_exp


# In[75]:


graph=pd.DataFrame({'Differences':Matches}, index=Probabilities)
graph.plot()
plt.show()


# In[22]:


Probabilities.append(P)
Matches.append(i)


# In[17]:


Probabilities=[]
Matches=[]


# In[ ]:


venn2([set(motifs_Hela), set(motifs_random)], set_labels = ('Background:HeLa', 'Background:Random'))
plt.title('Первичные мотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_single_Hela)):
    if (motifs_single_Hela[i] in motifs_single_random):
        print(motifs_single_Hela[i])


# In[ ]:





# In[91]:


ex=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(random)/table6_triple.csv', sep=',')
triple_random=ex.copy()
triple_random.rename(columns={'Unnamed: 0': 'Motifs'}, inplace=True)
del triple_random['P_ABC']
del triple_random['P_BC']
#double_random=double_random[:100]
#print(len(set(triple_random['Motifs'].values)))
new_column=[]
for elem in triple_random['Probability'].values:
    if elem!=math.inf:  
        i=round(-math.log(elem))
    else:
        i=0
    new_column.append(i)
triple_random['Score']=new_column
del triple_random['Probability']
cols = list(triple_random.columns)
a, b = cols.index('Score'), cols.index('Location')
cols[b], cols[a] = cols[a], cols[b]
triple_random = triple_random[cols]
i=[]
counter=[]
a=triple_random['Motifs'].values
for elem in a:
    if elem not in i:
        i.append(elem)
        counter.append(1)
    else:
        counter.append(0)
triple_random['Counter']=counter
triple_random=triple_random[triple_random['Counter']==1]
del triple_random['Counter']   
#triple_random=triple_random[:40]
triple_random.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/triple_random_all.csv')
triple_random
#motifs_triple_random=triple_random['Motifs'].values


# In[88]:


venn3([set(motifs_triple_HeLa), set(motifs_triple_proteom), set(motifs_triple_random)], set_labels = ('Background:HeLa', 'Background:Proteom', 'Background:Random'))
plt.title('Тройные мотивы для различных backgrounds (первые 10)\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_triple_HeLa)):
    if (motifs_triple_HeLa[i] in motifs_triple_proteom) and (motifs_triple_HeLa[i] in motifs_triple_random):
        print(motifs_triple_HeLa[i])


# In[244]:


len(single_random)


# In[96]:


#составляем единую таблицу для мотивов с бэкграундом HeLa
new_column=[]
single_random.sort_values('Probability',ascending=True,inplace=True)
for elem in single_random['Probability'].values:
    if elem!=0:
        i=round(-math.log(elem))
    else:
        i='P=0'
    new_column.append(i)
single_random['Score']=new_column
cols = list(single_random.columns)
a, b = cols.index('Score'), cols.index('Location')
cols[b], cols[a] = cols[a], cols[b]
single_random = single_random[cols]
del single_random['Probability']
single_random


# In[97]:


double_random=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/double_random_all.csv', sep=',')
del double_random['Unnamed: 0']
single_double_random=pd.merge(single_random,double_random,how='outer')
single_double_random


# In[98]:


single_double_triple_random=pd.merge(single_double_random,triple_random, how='outer')
single_double_triple_random.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/single_double_triple_random_all.csv')


# In[251]:


motifs_single_double_triple_HeLa=single_double_triple_HeLa['Motifs'].values
motifs_single_double_triple_proteom=single_double_triple_proteom['Motifs'].values
motifs_single_double_triple_random=single_double_triple_random['Motifs'].values


# In[252]:


venn3([set(motifs_single_double_triple_HeLa), set(motifs_single_double_triple_proteom), set(motifs_single_double_triple_random)], set_labels = ('Background:HeLa', 'Background:Proteom', 'Background:Random'))
plt.title('Mотивы для различных backgrounds\n')
print('Повторяющиеся мотивы')
for i in range(len(motifs_single_double_triple_HeLa)):
    if (motifs_single_double_triple_HeLa[i] in motifs_single_double_triple_proteom) and (motifs_single_double_triple_HeLa[i] in motifs_single_double_triple_random):
        print(motifs_single_double_triple_HeLa[i])


# In[173]:


#single_double_triple_HeLa
#ex=single_double_triple_HeLa.copy()
#ex.drop_duplicates()
print(len(set(motifs_single_double_triple_HeLa)))


# In[99]:


single_double_triple_HeLa_1=single_double_triple_HeLa.copy()
single_double_triple_HeLa_1=single_double_triple_HeLa_1.set_index('Motifs')
single_double_triple_proteom_1=single_double_triple_proteom.copy()
single_double_triple_proteom_1=single_double_triple_proteom_1.set_index('Motifs')
single_double_triple_random_1=single_double_triple_random.copy()
single_double_triple_random_1=single_double_triple_random_1.set_index('Motifs')


# In[254]:


result=pd.merge(single_double_triple_HeLa_1,single_double_triple_proteom_1, how='inner', on=index)
result


# In[102]:


result=pd.merge(single_double_triple_HeLa_1,single_double_triple_proteom_1,how='outer', left_index=True, right_index=True)
del result['Occurrences_y']
del result['Location_x']
del result['Location_y']
result.rename(columns={'Score_x': 'Score_HeLa','Score_y': 'Score_proteom','Occurrences_x': 'Occurrences','Backgrounds_x':'Backgrounds_HeLa','Backgrounds_y':'Backgrounds_proteom'}, inplace=True)
cols = list(result.columns)
a, b = cols.index('Score_HeLa'), cols.index('Occurrences')
cols[b], cols[a] = cols[a], cols[b]
result


# In[105]:


result_final=pd.merge(single_double_triple_random_1,result, how='outer',left_index=True, right_index=True)
del result_final['Occurrences_y']
del result_final['Location']
result_final.rename(columns={'Score': 'Score_random','Backgrounds': 'Backgrounds_random','Occurrences_x': 'Occurrences'}, inplace=True)
result_final.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/results_final_all.csv')


# In[106]:


q=[]
for elem in result_final['Score_random'].values:
    if elem=='P=0':
        q.append(10000)
    else:
        q.append(elem)
result_final['Score_random']=q
result_final


# In[107]:


result_final.sort_values('Score_random',ascending=False,inplace=True)
q=[]
for elem in result_final['Score_random'].values:
    if elem==10000:
        q.append('P=0')
    else:
        q.append(elem)
result_final['Score_random']=q
result_final


# In[195]:


i=[]
counter=[]
a=result_final.index
for elem in result_final.index:
    if elem not in i:
        i.append(elem)
        counter.append(1)
    else:
        counter.append(0)
        
print(len(i))        


# In[108]:


#result_final['Counter']=counter
#result_final=result_final[result_final['Counter']==1]
#del result_final['Counter']
result_final.to_csv('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/FINAL_new_all.csv')
result_final.to_excel("/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/FINAL_new_all.xlsx")


# In[260]:


cols = list(result_final.columns)
cols[1],cols[3], cols[5], cols[0], cols[4], cols[6], cols[2]=cols[0], cols[1], cols[2],cols[3], cols[4], cols[5], cols[6]
result_final


# In[264]:


result_final.to_excel("/home/vikalin/Article/example/PHOSPHORYLATION/S_6_backgrounds_results/FINAL_new.xlsx")


# In[262]:


result_final.drop_duplicates()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[74]:


print(P[9][:])


# In[ ]:





# In[76]:


x=sum(1 for interval in intervals if (interval[7]=='P'))
print('x:',x)


# In[71]:


x=sum(1 for interval in intervals if (( (interval[0]=='R') 
                                       and (interval[7]=='P'))))

print(x)

u=sum(1 for interval in background if (( (interval[0]=='R')                                        
                                       and (interval[7]=='P'))))

print(u)


# In[95]:


example=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/schwartz_triple_last.csv', sep=',')
single=example.copy()
-math.log(single.iloc[90]['Probability'])


# In[99]:



proteins=pd.read_csv('/home/vikalin/20151001_03_Qp1_HeLa_1ug.mzML.demix_proteins.csv',sep='\t')
hela_protein=proteins['dbname'].values
hela_protein


# In[101]:


proteins=pd.read_csv('/home/vikalin/20151001_03_Qp1_HeLa_1ug.mzML.demix_proteins.csv',sep='\t')
hela_protein=proteins['dbname'].values
hela=dict()
for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):    
    #ищем его номер
    i=description.find('|')
    k=description.rfind('|')
    s=description[i+1:k]
    if s in hela_protein:
        hela[s]=sequence
        print('!')


# In[102]:


print(len(hela.keys()))


# In[119]:


score=[]
for i in range(len(single['Probability'].values )):
    if single['Probability'].values[i]!=0:
        score.append(-math.log(single['Probability'].values[i]))
    else:
        score.append(None)      
single['Score']=score 
single=single.dropna() 
single


# In[135]:


triple_motifs=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/schwartz_triple.csv', sep=',')
triple_motifs.sort_values('Probability',ascending=True,inplace=True)
triple_motifs.where(triple_motifs['Probability']<10**(-16)).dropna()[::2]


# In[133]:


double_motifs=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/schwartz.csv', sep=',')
double_motifs.where(double_motifs['Probability']<10**(-20)).dropna()


# In[121]:


k=0
for interval in intervals:
    if interval[1]=='R' and interval[3]=='R' and interval[4]=='S':
        k+=1
print(k)
print(len(intervals))


# In[122]:


k=0
for elem in background:
    if elem[1]=='R' and elem[3]=='R' and elem[4]=='S':
        k+=1
print(k)
print(len(background))


# In[115]:


P[18][7]


# In[99]:


print(P)


# In[54]:


type(result['Список мотивов'].values)


# In[132]:


table = pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/schwartz.csv', sep=',')
table=table.where(table['First acid']=='E')
table=table.dropna()
table[:40]


# In[151]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
double_motifs=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/schwartz.csv', sep=',')
primary_motifs=pd.read_csv('/home/vikalin/Article/example/PHOSPHORYLATION/primary.csv', sep=',')
triple_motifs=triple_motifs_creator(acids,double_motifs,primary_motifs,background,intervals,'S',6,'/home/vikalin/Article/example/PHOSPHORYLATION/schwartz_triple.csv')


# In[108]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
occurrences[15][4]


# In[ ]:





# In[ ]:





# In[ ]:


file=open('/home/vikalin/knowledgebase/PTM/PTM.txt','r')
for line in file:
    l=line.split('\t')
    print(l[0],l[2])
    break


# In[ ]:


path_FASTA='/home/vikalin/knowledgebase/proteomes/HUMAN.fasta.gz'
a=[]
for description, sequence in fasta.read(gzip.open(path_FASTA,'rt')):
    a.append(description)
print(len(a))    


# In[ ]:


'knowledgebase/proteomes/HUMAN.fasta.gz'


# In[ ]:


interval_length=6
modification='PHOSPHORYLATION'
modification_site='S'
acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
intervals=mod_intervals_DB('/home/vikalin/Article/identification.txt','/home/vikalin/Article/HUMAN.fasta.gz',Peptides,interval_length,modification_site)


# In[ ]:


background=background_maker(modification_site,interval_length,'/home/vikalin/knowledgebase/PTM/HUMAN_proteinDB.fasta.gz','/home/vikalin/knowledgebase/PTM/human_phosphorylation.txt',modification)


# In[ ]:


P=P_counter(intervals,interval_length,modification_site,acids,background)


# In[ ]:


heatmap_visualization(P,acids,interval_length,modification_site,'Исходные вероятности для '+modification+' '+modification)


# In[ ]:


occurrences=occurrences_counter(intervals,interval_length,acids,modification_site)


# In[ ]:


heatmap_visualization(occurrences,acids,interval_length,modification_site,'Occurrences для '+modification+' '+modification_site)


# In[ ]:


volcano_plot(acids,interval_length,occurrences,P,modification_site,modification)


# In[ ]:





# In[61]:


table = pd.read_csv('/home/vikalin/Article/example/table15.csv', sep=',')
table['Background'].values[10]


# In[7]:


ex=open('/home/vikalin/intervals.txt','r')
intervals=[]
for line in ex:
    line=line.replace('\n','')
    intervals.append(line)


# In[11]:


i=0
for elem in intervals:
    if elem[7]=='P':
        i+=1
print(i,len(intervals),(i/len(intervals))*100)
        


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[2]:


ex=open('/home/vikalin/Article/example/PHOSPHORYLATION/S_6_background:FASTA(4 буквы+хи квадрат)/occurrences.txt','r')
chi2_results=[]
for i in range(0,21):
    a=[0]*(6*2+1)
    chi2_results.append(a)
k=0
for line in ex:
    s=line.split(' ')
    for i in range(len(s)):
        if s[i]!='\n':
            chi2_results[k][i]=int(s[i])
    k+=1
print(chi2_results)    


# In[12]:


acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']

farmers = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]




fig, ax = plt.subplots(figsize=(15, 15))
im = ax.imshow(chi2_results)

# We want to show all ticks...
ax.set_xticks(np.arange(len(farmers)))
ax.set_yticks(np.arange(len(acids)))
# ... and label them with the respective list entries
ax.set_xticklabels(farmers)
ax.set_yticklabels(acids)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels())

# Loop over data dimensions and create text annotations.
for i in range(len(acids)):
    for j in range(len(farmers)):
        text = ax.text(j, i, chi2_results[i][j],
                       ha="center", va="center", color="w",fontsize=17)

ax.set_title("Встречаемость аминокислот")
plt.show()


# In[ ]:




