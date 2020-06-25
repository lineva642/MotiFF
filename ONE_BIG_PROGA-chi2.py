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

path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
background_type='FASTA'
# print('Название исследуемого набора данных:')
name_sample='example'
# print('Название организма:')
organism='Human'
# print('Введите длину интервала:')
interval_length=6
# print('Введите модификацию:')
modification='PHOSPHORYLATION'
# print('Введите сайт модификации:')
modification_site='S'


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
print('Peptides',Peptides)



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
def mod_intervals_DB_experimental_2(path_FASTA,Peptides,interval_length,modification_site,path_sample):    
    #mod_intervals=interval_maker_experimental(idPeptides['Protein'].values,idPeptides['Peptide'].values,idPeptides['Mod_Peptide'].values,interval_length,modification_site)
    mod_intervals=interval_maker_experimental(Peptides,interval_length,modification_site)
    print('Intervals DB is ready!')
    return mod_intervals





# In[93]:


# In[94]:


#набор функции для построения базы немодифицированных интервалов или backgroundов


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
def FASTA_DB_creation_random(path_FASTA):
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
def FASTA_DB_creation_HeLa(path_FASTA):
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


def background_maker(modification_site,interval_length,path_FASTA,modification,intervals):
    print('Making background DB')
    background=background_array_FASTA(path_FASTA,modification_site,interval_length,intervals)
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
    saving_table(path_results,primary_motif,interval_length,'primary')
    print('Primary motifs are ready!')
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
    
    saving_table(path_results,second_motifs_selection_p_copy,interval_length,'double')
    print('Double motifs are ready!')
    
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
    print('Triple motifs are ready!')
    
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
    
    print('Quadro motifs are ready!')
    
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



# In[118]:


#организация хранения результатов
def saving(name_sample,interval_length,modification,modification_site): 
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


def output_experimental(name_sample,Peptides,interval_length,modification_site,modification,path_FASTA):
    #path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
    #path_identification='/home/vikalin/Article/identification.txt'
    path_sample,path_modification,path_results=saving(name_sample,interval_length,modification,modification_site)
    acids=['Y','W','F','M','L','I','V','A','C','P','G','H','R','K','T','S','Q','N','E','D',' ']
    intervals=mod_intervals_DB_experimental(path_FASTA,Peptides,interval_length,modification_site,path_modification)
    print('Длина интервалов в exp',len(intervals))
    background=background_maker(modification_site,interval_length,path_FASTA,modification,intervals)
    print('Длина интервалов в back',len(background))
    occurrences=n_calculator(acids,intervals,interval_length,path_results)
    background_n=background_n_matrix(acids,interval_length,background,path_results)
    expected_FASTA,chi2_results,chi2_selection=p_value(background_n,occurrences,interval_length,modification_site,acids,path_results)
    single,double,triple,quadruple=motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,path_results)
    print(single,double,triple,quadruple)
#     single,double,triple,quadruple=motifs(acids,chi2_selection,interval_length,modification_site,background,intervals,path_results)
#     print(single,double,triple,quadruple)
    
    return chi2_results,chi2_selection,intervals,background

# In[122]:


#поверка на наборе данных


# In[123]:


# #загружаем таблицу модифицированных пептидов
# if not os.path.isfile('/home/vikalin/Article/dataset.xls'):
#     print ('Downloading the dataset')
#     urlretrieve(
#         'https://www.pnas.org/content/pnas/suppl/2004/08/06/0404720101.DC1/04720Table3.xls',
#         '/home/vikalin/Article/dataset.xls')
#     print ('Done!')
# else:
#     print('XLS file for dataset already exists!')


# # In[124]:


# #улучшение проги
# xl = pd.ExcelFile('/home/vikalin/Article/dataset.xls')
# df1 = xl.parse('Sheet1')
# df1.columns=df1.loc[0]
# df=df1.drop(0)
# Peptides=((df.loc[df['Amb?'] == 'UNIQUE']).reset_index()).loc[:,df.columns.intersection(['Peptide']) ]
# indexes=[elem.replace('*','') for elem in Peptides['Peptide']]
# Peptides['Mod_Peptide']=indexes  
# Peptides.columns = ['Mod_Peptide','Peptide']
# Peptides['Protein']=None
# ind=Peptides['Peptide'].values
# Peptides['index']=ind
# Peptides=Peptides.set_index('index')
# Peptides



# In[125]:


# path_FASTA='/home/vikalin/Article/HUMAN.fasta.gz'
# background_type='FASTA'
# print('Название исследуемого набора данных:')
# name_sample=str(input())
# print('Название организма:')
# organism=str(input())
# print('Введите длину интервала:')
# interval_length=int(input())
# print('Введите модификацию:')
# modification=str(input())
# print('Введите сайт модификации:')
# modification_site=str(input())
chi2_results,chi2_selection,intervals,background=output_experimental(name_sample,Peptides,interval_length,modification_site,modification,path_FASTA)