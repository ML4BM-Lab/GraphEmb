import numpy as np
import pandas as pd
import random
import os
from rdkit import Chem
import requests
random.seed(1)

def getamino_KEGG(protein):
    r = requests.get(f'http://rest.kegg.jp/get/{protein}/aaseq')
    aminoseq = ''.join(r.text.split('\n')[1:])
    print(aminoseq)

    if aminoseq =='':
        return None
    return aminoseq

def get_drug_KEGG(drug):
    r = requests.get(f'http://rest.kegg.jp/get/{drug}/mol')
    mol = Chem.MolFromMolBlock(r.text)
    smiles = Chem.MolToSmiles(mol)
    print(smiles)
    return smiles

def get_dataset(name):
    #Read csv to dataframe for positive pairs
    colnames = ['Drug ID', 'Gene']
    data_path =os.getcwd() + '/../../DB/Data/Yamanashi_et_al_GoldStandard/' + name.upper()+'/interactions/' + name+'_admat_dgc_mat_2_line.txt'
    df = pd.read_csv(data_path, header = None, names = colnames, index_col=False, sep='\t')

    df['Label'] = 1
    #example: hsa10 -> hsa:10
    df['Gene'] = df['Gene'].map(lambda x: x[:3] + ":" + x[3:])

    #Obtain the unique drugs and genes
    drugs = df['Drug ID'].unique()
    genes = df['Gene'].unique()

    #print(drugs)

    #Get the smiles and sequences of drugs and proteins
    fname = 'drugs_smiles_Yamanishi_' + name + '.txt'
    path_smiles = os.getcwd() + '/../Data/Yamanishi_et_al_GoldStandard/' + fname
    if not os.path.exists(path_smiles):
        smiles = np.vectorize(get_drug_KEGG)(drugs)
        df_smiles = pd.DataFrame()
        df_smiles['Drug ID'] = drugs
        df_smiles['SMILES'] = smiles
        df_smiles.to_csv(path_smiles, header=None, index = None, sep = ' ')
        
    drugs_smiles = pd.read_csv(path_smiles, delimiter=" ", header = None)
    drugs_smiles = dict(zip(drugs_smiles[0], drugs_smiles[1]))

    fname = 'genes_targetsequences_Yamanishi_' + name + '.txt'
    path_targetsequences = os.getcwd() + '/../Data/Yamanishi_et_al_GoldStandard/' + fname
    if not os.path.exists(path_targetsequences):
        targetsequences = np.vectorize(getamino_KEGG)(genes)
        df_genes = pd.DataFrame()
        df_genes['Gene'] = genes
        df_genes['Target Sequence'] = targetsequences
        df_genes.to_csv(path_targetsequences, header=None, index = None, sep = ' ')
    
    genes_targetsequences = pd.read_csv(path_targetsequences, delimiter=" ", header = None)
    genes_targetsequences = dict(zip(genes_targetsequences[0], genes_targetsequences[1]))   

    #print(smiles)

    #Create dictionaries
    #drugs_smiles = {d:s for d, s in zip(drugs, smiles)}
    #genes_targetsequences = {g:t for g,t in zip(genes, targetsequences)}
    #Add columns for corresponding SMILES and Target Sequence of each pair
    df['SMILES'] = df['Drug ID'].map(drugs_smiles)
    df['Target Sequence'] = df['Gene'].map(genes_targetsequences)
    print(len(df.index))
    df = df.replace(to_replace='None', value=np.nan).dropna()
    print(len(df.index))

    #Get data for splits
    columns = df[['Drug ID','Gene', 'Label']]
    df_splits = columns.copy()
    output_path=os.getcwd() + '/../Data/Yamanishi_et_al_GoldStandard/Yamanishi_' + name + '_pairs.txt'
    df_splits.to_csv(output_path)


    #Randomnly choose unseen pairs
    seenpairs = set(df['Drug ID'] + df['Gene'])
    unseenpairs = set()
    n = len(seenpairs)
    while n > 0:
        drug_sample = random.choice(drugs)
        gene_sample = random.choice(genes)
        
        if not genes_targetsequences[gene_sample]=='None':
            result = drug_sample + gene_sample
            if (result not in seenpairs and result not in unseenpairs):
                sample = pd.DataFrame(data = {'Drug ID':[drug_sample], 'Gene': [gene_sample], 'Label':[0]})
                df = pd.concat([df, sample], ignore_index = True)
                unseenpairs.add(result)
                n-=1
        

    #Add columns for corresponding SMILES and Target Sequence of each pair
    df['SMILES'] = df['Drug ID'].map(drugs_smiles)
    df['Target Sequence'] = df['Gene'].map(genes_targetsequences)
   


    #Save it as a csv file
    output_path = os.getcwd() + '/../Data/Yamanishi_et_al_GoldStandard/Yamanishi_' + name.upper() + ".csv"
    df.to_csv(output_path)

    print("Yamanishi " + name + " done")


get_dataset('e')
get_dataset('gpcr')
get_dataset('ic')
get_dataset('nr')