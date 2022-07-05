import pandas as pd
import numpy as np
import os
import sys
import functools as ft
from tqdm import tqdm
from numpy import savetxt

def format_matrices(dataset,subdataset):

    #define path
    fpath = os.path.join(os.getcwd(),'NeoDTI/Data',dataset,subdataset)
    
    #print(f"path {fpath}")

    def retrieve_matrices():
        
        prot_f = []
        drug_f = []

        #categorize in drug/proteins files
        for f in os.listdir(fpath):
            if 'mat_prot' in f or 'Proteins' in f:
                prot_f.append(f)
            elif 'mat_drug' in f or 'Drugs' in f:
                drug_f.append(f)

        #return the sorted array
        return sorted(prot_f), sorted(drug_f)

    def retrieve_dti():

        dtipath = [f for f in os.listdir(fpath) if '2line' in f][0]

        #read the f
        dtif = pd.read_csv(os.path.join(fpath,dtipath),sep='\t',header=None)
        dtif.columns = ['Drugs','Targets']

        return dtif.Drugs.tolist(), dtif.Targets.tolist()

    def filter_matrices_proteins_symmat(mat_l,savepath,dtil=None):

        names = []

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)

            try:
                matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            except:
                matnames = mat.columns.to_list()

            print(f'{matnames[0:5]} for file {f} (targets)')
            names.append(matnames)

        #now select the intersection
        inames = sorted(ft.reduce(lambda a,b: set(a).intersection(set(b)),names))

        if dtil:
            inames = sorted(set(inames).intersection(set(dtil)))

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)

            try:
                matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            except:
                matnames = mat.columns.to_list()
            
            mat.columns = matnames
            mat.index = matnames
            #generate the final mat and save
            fmat = mat.loc[inames,inames]
            print(f"shape of formatted target matrix-> {fmat.shape}")
            fmat.to_csv(os.path.join(savepath,'Form_symmat_preSNF_'+f),sep='\t')

    def filter_matrices_drugs_symmat(mat_l,savepath,dtil=None):

        def sortedIntersect(v1,v2):
            sortInt = []
            for itemV1 in v1:
                if itemV1 in v2:
                    sortInt.append(itemV1)
            return sortInt

        names = []

        #get the MorganFingerprint file
        for f in os.listdir(fpath):
            if 'MorganFingerprint' in f:
                Morganf = f

        

        #read the morganfingerprint file
        morganMat = pd.read_csv(os.path.join(fpath,Morganf),sep='\t',index_col=0)
        morganIndex = mat.index.tolist()
        refMat = pd.read_csv(os.path.join(fpath,'drug.txt'),header=None).values.flatten().tolist()
        sortedIntersect = [drug for drug in refMat if drug in morganIndex]

        #now select the intersection
        inames = sorted(ft.reduce(lambda a,b: set(a).intersection(set(b)),names))

        if dtil:
            print(f"inames type {type(inames[0])}, dtil type {type(dtil[0])}")
            inames = sorted(set(inames).intersection(set(dtil)))

        for f in mat_l:

            if 'aers' in f and dataset == 'BindingDB':
                mat = pd.read_csv(os.path.join(fpath,f),sep=',',index_col=0)
            else:
                if 'Rchem' not in f:
                    mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
                else:
                    mat = pd.read_csv(os.path.join(fpath,f),sep=' ',index_col=0)

            try:
                matnames = [drug.replace(':','') for drug in mat.index.to_list()]
            except:
                matnames = mat.index.to_list()

            mat.index = matnames
            mat.columns = matnames
            
            #intersect names
            fmat = mat.loc[inames,inames]

            print(f"shape of formatted drug matrix-> {fmat.shape}")

            fmat.to_csv(os.path.join(savepath,'Form_symmat_preSNF_'+f),sep='\t')

    
    prot_f, drug_f = retrieve_matrices()
    #dtidrug, dtiprot = retrieve_dti()

    #print(f'Matrices for proteins {prot_f}')
    #print(f'Matrices for drugs {drug_f}')

    folder_path = os.path.join(fpath,'Formatted/')
    #print(f"folder path {folder_path}")

    if not os.path.isdir(folder_path):
        #create the folder
        os.mkdir(os.path.join(fpath,'Formatted'))

    #------------------- Format matrices ------------------ #

    print("\ngenerating symmat")
    filter_matrices_proteins_symmat(prot_f,savepath = folder_path)
    filter_matrices_drugs_symmat(drug_f,savepath = folder_path)
    

#get sys
try:
    arg1 = sys.argv[1]
    try: 
        arg2 = sys.argv[2]
    except:
        arg2 = ''
except:
    arg1, arg2 = '','' 

dataset, subdataset = arg1, arg2


format_matrices(dataset,subdataset)
