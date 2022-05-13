from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import os
from tqdm import tqdm
import numpy as np
import sys
import pandas as pd

def MorganDiceSimilarity(mol1,mol2,moldd):
    return DataStructs.DiceSimilarity(moldd[mol1],moldd[mol2])

def MorganAdjMat(molpath,savepath):

    #check if file exists
    if os.path.exists(savepath+'MorganFingerprint.tsv'):
        print('File already generated!')
        return 

    moldd = {}

    for fname in tqdm(os.listdir(molpath), 'Building Morgan Dict'):
        try:
            molf = Chem.MolFromMolFile(os.path.join(molpath,fname))
            moldd[fname[:-4]] = AllChem.GetMorganFingerprint(molf,2)
        except:
            print(f'Could not compute Morgan Fingerprint for Drug {fname}')

    GenerateMatrix(moldd,savepath)

def GenerateMatrix(moldd,savepath):

    drugnames = sorted(moldd.keys())
    drugnames_L = len(drugnames)
    adjmat = np.zeros(shape=(drugnames_L,drugnames_L))

    for i in tqdm(range(drugnames_L), 'Building Dice Similarity'):
        for j in range(i,drugnames_L):
            adjmat[i,j] = MorganDiceSimilarity(drugnames[i],drugnames[j],moldd)

    #flip and fill the matrix (its symmetric)
    utriindex = np.triu_indices(adjmat.shape[0],k=1)
    adjmat.T[utriindex] = adjmat[utriindex]

    #to df and save
    adjmat_df = pd.DataFrame(adjmat,index = drugnames, columns = drugnames)
    adjmat_df.to_csv(savepath+'MorganFingerprint.tsv', sep='\t')

molpath, savepath = sys.argv[1], sys.argv[2]
MorganAdjMat(molpath, savepath)