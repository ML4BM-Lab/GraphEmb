import numpy as np
import pandas as pd
import random
import os
from pubchempy import Compound
import requests
random.seed(1)

def getamino_KEGG(protein):
    r = requests.get(f'http://rest.kegg.jp/get/{protein}/aaseq')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

def get_drug_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles



def get_dataset(type):
    #Read csv to dataframe for positive pairs
    colnames = ['Drug ID Kegg', 'Gene']
    data_path =os.getcwd() + '/../../../Data/Yamanashi_et_al_GoldStandard/' + type.upper()+'/interactions/' + type+'_admat_dgc_mat_2_line.txt'
    df = pd.read_csv(data_path, header = None, names = colnames, index_col=False, sep='\t')

    colnames = ['pubchem','dr']
    data_path = os.getcwd()+ '/../Data/Yamanashi_et_al_GoldStandard/pubchem.txt'
    pubchem_dr = pd.read_csv(data_path, header = None, names = colnames, sep='\t')

    pubchem_dr['pubchem'] = pubchem_dr['pubchem'].map(lambda x: x.lstrip('pubchem:'))
    pubchem_dr['dr'] = pubchem_dr['dr'].map(lambda x: x.lstrip('dr:'))

    pubchem_dr_dict = dict(zip(pubchem_dr['dr'], pubchem_dr['pubchem']))

    df['Drug ID'] = df['Drug ID Kegg'].map(pubchem_dr_dict)

    df = df.drop(columns = ['Drug ID Kegg'])

    df['Label'] = 1

    df['Gene'] = df['Gene'].map(lambda x: x[:3] + ":" + x[3:])
    print(df)

    #Obtain the unique drugs and genes
    drugs = df['Drug ID'].unique()
    genes = df['Gene'].unique()

    #print(drugs)

    #Get the smiles and sequences of drugs and proteins
    fname = 'smiles_Yamanashi_' + type + '.txt'
    path_smiles = os.getcwd() + '/../Data/Yamanashi_et_al_GoldStandard/' + fname
    
    if not os.path.exists(path_smiles):
        smiles = np.vectorize(get_drug_pubchem)(drugs)
        np.savetxt(path_smiles, smiles, fmt="%s")
    else:
        smiles = np.genfromtxt(path_smiles, dtype = 'str')
    
    fname = 'targetsequences_Yamanashi_' + type + '.txt'
    path_targetsequences = os.getcwd() + '/../Data/Yamanashi_et_al_GoldStandard/' + fname
    
    if not os.path.exists(path_targetsequences):
        targetsequences = np.vectorize(getamino_KEGG)(genes)
        np.savetxt(path_targetsequences, targetsequences, fmt="%s")
    
    else:
        targetsequences = np.genfromtxt(path_targetsequences, dtype = 'str')

    #print(smiles)

    #Create dictionaries
    drugs_smiles = {d:s for d, s in zip(drugs, smiles)}
    genes_targetsequences = {g:t for g,t in zip(genes, targetsequences)}

    #Randomnly choose unseen pairs
    seenpairs = set(df['Drug ID'] + df['Gene'])
    unseenpairs = set()
    n = len(seenpairs)
    while n > 0:
        drug_sample = random.choice(drugs)
        gene_sample = random.choice(genes)
        
        result = drug_sample + gene_sample
        if (result not in seenpairs and result not in unseenpairs):
            sample = pd.DataFrame(data = {'Drug ID':[drug_sample], 'Gene': [gene_sample], 'Label':[0]})
            df = pd.concat([df, sample], ignore_index = True)
            unseenpairs.add(result)
            n-=1
        

    #Add columns for corresponding SMILES and Target Sequence of each pair
    df['SMILES'] = df['Drug ID'].map(drugs_smiles)
    df['Target Sequence'] = df['Gene'].map(genes_targetsequences)    
    print(df)

    print(df.columns)
    cols = ['Drug ID', 'Gene', 'SMILES', 'Target Sequence', 'Label']
    df = df[cols]

    #Save it as a csv file
    output_path = os.getcwd() + '/../Data/Yamanashi_et_al_GoldStandard/Yamanashi_' + type.upper() + ".txt"
    df.to_csv(output_path, header=None, index = None, sep = ' ')


get_dataset('e')
get_dataset('gpcr')
get_dataset('ic')
get_dataset('nr')