import pandas as pd
import numpy as np
import os
import sys
import functools as ft
from tqdm import tqdm
from numpy import savetxt
from collections import defaultdict as dd

def format_matrices(dataset,subdataset):

    #define path
    fpath = os.path.join(os.getcwd(),'NeoDTI/Data',dataset,subdataset)

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

    def generate_dti_mat(drug_shape, prot_shape):

        def dict_from_dti(dti):
            dtidd= dd(int)
            for pair in dti:
                dtidd['||'.join(pair)] = 1
            return dtidd

        drugs = pd.read_csv(os.path.join(fpath,'drug.txt'), header=None).values.flatten().tolist()
        prots = pd.read_csv(os.path.join(fpath,'protein.txt'), header=None).values.flatten().tolist()

        print(f"Total drugs -> {len(drugs)}, drugs used {drug_shape}")
        print(f"Total proteins -> {len(prots)}, proteins used {prot_shape}")

        if len(drugs) == drug_shape and len(prots) == prot_shape:


            print("Consistency accomplished!")
            #read the f
            dtipath = [f for f in os.listdir(fpath) if 'DTI' in f][0]
            dtif = pd.read_csv(os.path.join(fpath,dtipath),sep='\t').values
            dtidd = dict_from_dti(dtif)

            #define matrix
            mat_drug_mat = np.zeros(shape=(drug_shape,prot_shape),dtype=int)

            for i,drug in enumerate(drugs):
                for j,prot in enumerate(prots):

                    mat_drug_mat[i,j] = dtidd['||'.join([drug,prot])]

            np.savetxt(os.path.join(fpath,'mat_drug_protein.txt'),mat_drug_mat,delimiter='\t',fmt='%i')

        else:
            print(f"Inconsistency!\n Drugs: {drugs} \n Proteins: {prots}")

    def filter_matrices_proteins_symmat(mat_l):

        shapes_v = []

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',header=None)
            shapes_v.append(mat.shape[0])

        if len(set(shapes_v)) == 1:
            return shapes_v[0]
        else:
            print("Inconsistency of samples")
            
    def filter_matrices_drugs_symmat(mat_l):

        shapes_v = []

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',header=None)
            shapes_v.append(mat.shape[0])

        if len(set(shapes_v)) == 1:
            return shapes_v[0]
        else:
            print("Inconsistency of samples")

    folder_path = os.path.join(fpath,'Formatted/')

    if not os.path.isdir(folder_path):
        #create the folder
        os.mkdir(os.path.join(fpath,'Formatted'))

    #------------------- Format matrices ------------------ #

    prot_f, drug_f = retrieve_matrices()
    drug_shape = filter_matrices_drugs_symmat(drug_f)
    prot_shape = filter_matrices_proteins_symmat(prot_f)
    generate_dti_mat(drug_shape, prot_shape)
    

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
