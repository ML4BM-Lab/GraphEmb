import numpy as np
import pandas as pd
import random
import os
from pubchempy import Compound
import requests
random.seed(1)


def get_drug_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

def getamino_uniprot(protein):
    r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

###----------------------------------Yamanishi-------------------------------------###
#Read csv to dataframe for positive pairs
df = pd.read_csv(os.getcwd() + '/interactions/All_DB_Did_Tid_HumanSingleProteins_FDA_D_withSMI_T_NoFalseAA.txt', index_col=False, sep='\t')
#df = df.rename(columns ={'#Drug':'DrugBank ID'})
df['Label'] = [1]*len(df.index)
df.columns = ['DrugBank ID' , 'Gene', 'Label']

#Obtain the unique drugs and genes
drugs = df['DrugBank ID'].unique()
genes = df['Gene'].unique()

print(drugs)

#Get the smiles and sequences of drugs and proteins
smiles = np.vectorize(getdrug_drugbank)(drugs)
targetsequences = np.vectorize(getamino_uniprot)(genes)

print(smiles)

#Create dictionaries
drugs_smiles = {d:s for d, s in zip(drugs, smiles)}
genes_targetsequences = {g:t for g,t in zip(genes, targetsequences)}

#Randomnly choose unseen pairs
seenpairs = set(df['DrugBank ID'] + df['Gene'])
unseenpairs = set()
n = len(seenpairs)
while n > 0:
    drug_sample = random.choice(drugs)
    gene_sample = random.choice(genes)
    
    result = drug_sample + gene_sample
    if (result not in seenpairs and result not in unseenpairs):
        sample = pd.DataFrame(data = {'DrugBank ID':[drug_sample], 'Gene': [gene_sample], 'Label':[0]})
        df = pd.concat([df, sample], ignore_index = True)
        unseenpairs.add(result)
        n-=1
    

#Add columns for corresponding SMILES and Target Sequence of each pair
df['SMILES'] = df['DrugBank ID'].map(drugs_smiles)
df['Target Sequence'] = df['Gene'].map(genes_targetsequences)    

#Save it as a csv file
df.to_csv('DrugBankFDA.csv')