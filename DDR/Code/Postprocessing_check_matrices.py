import pandas as pd
import numpy as np
import os

def check_diag(m):
    print(f'Diag -> {np.diag(m)[0:5]}')
    print(f'Checking diag to 1 -> {all(np.diag(m) == np.ones(m.shape[0]))}')

def check_simmetry(m):

    #get indices
    utriindex = np.triu_indices(m.shape[0],k=1)

    #get the values
    uppertri = np.array(m)[utriindex]
    lowertri = np.array(m.T)[utriindex]

    print(f'Upper tri -> {uppertri[0:5]}')
    print(f'Lower tri -> {lowertri[0:5]}')

    #check if 
    print(f'Checking simmetry -> {all(uppertri == lowertri)}')

def check_mat(path):
    #Lambda
    print('Loading mat')
    lambda_mat = pd.read_csv(path,index_col=0,sep='\t')
    #check diag to 1
    check_diag(lambda_mat)
    #check simmetry
    check_simmetry(lambda_mat)

def check_properties():

    matrices = [f'drug_Rchemcpp_{mat}' for mat in ['lambda','marginalized','minmaxTanimoto','spectrum','tanimoto']] +\
               ['Prot_BioGrid_SHD','Prot_GO_PPI']

    for DB in ['BindingDB','BIOSNAP','Davis_et_al','DrugBank','Yamanashi_et_al_GoldStandard']:
        print(f'\n\nDataset: {DB}')
        if DB == 'Yamanashi_et_al_GoldStandard':
            for DF in ['E','GPCR','IC','NR']:
                print(f'\nDataframe: {DF}')
                for mat in matrices:
                        print(f'\nMat: {mat}')
                        check_mat(path=os.getcwd()+f'/DDR/Data/{DB}/{DF}/{DF}_{mat}.tsv')
        else:
            for mat in matrices:
                    print(f'\nMat: {mat}')
                    check_mat(path=os.getcwd()+f'/DDR/Data/{DB}/{DB}_{mat}.tsv')

print('Checking properties!')
check_properties()