import pandas as pd
import numpy as np
import os
import sys
import functools as ft

def format_matrices(dataset,subdataset):

    #define path
    fpath = os.path.join(os.getcwd(),'DDR/Data',dataset,subdataset)

    def retrieve_matrices():
        
        prot_f = []
        drug_f = []

        #categorize in drug/proteins files
        for f in os.listdir(fpath):
            if 'prot' in f:
                prot_f.append(f)
            elif 'drug' in f:
                drug_f.append(f)

        #return the sorted array
        return sorted(prot_f), sorted(drug_f)

    def filter_matrices_proteins(mat_l):

        names = []

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
            matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            names.append(matnames)

        #now select the intersection
        inames = sorted(ft.reduce(lambda a,b: set(a).intersection(set(b)),names))

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
            matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            mat.columns = matnames
            mat.index = matnames
            #generate the final mat and save
            fmat = mat.loc[inames,inames]
            fmat.to_csv(os.path.join(fpath,'Formatted','Form_'+f),sep='\t')

    def filter_matrices_drugs(mat_l):

        names = []

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
            
            if 'aers' in f:
                matnames = mat.index.to_list()
                names.append(matnames)

        #now select the intersection
        inames = sorted(ft.reduce(lambda a,b: set(a).intersection(set(b)),names))

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
            matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            mat.columns = matnames
            mat.index = matnames
            #generate the final mat and save
            fmat = mat.loc[inames,inames]
            fmat.to_csv(os.path.join(fpath,'Formatted','Form_'+f),sep='\t')


    prot_f, drug_f = retrieve_matrices()

    #create the folder
    os.mkdir(os.path.join(fpath,'Formatted'))

    #save the new formatted matrices
    filter_matrices_proteins(prot_f)


    
