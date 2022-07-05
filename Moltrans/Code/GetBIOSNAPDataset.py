import numpy as np
import pandas as pd
import random
import requests
import os
from getSmilesDrugBank import *
random.seed(1)

###-----------------------------------------------BIOSNAP-------------------------------------------##
#The drugs can be obtained from the link 
#https://go.drugbank.com/structures/small_molecule_drugs/{DRUG}.smiles
#Example -> Drug: DB00755 -> https://go.drugbank.com/structures/small_molecule_drugs/DB00755.smiles

#The proteins can be obtained from the uniprot link
#https://www.uniprot.org/uniprot/{PROTEIN}.fasta
#Example -> Protein: Q16539 -> https://www.uniprot.org/uniprot/Q16539.fasta

def getdrug_drugbank(drug):
    r = requests.get(f'https://go.drugbank.com/structures/small_molecule_drugs/{drug}.smiles')
    return r.text

def getamino_uniprot(protein):
    r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    if aminoseq == '':
        return (protein, None)
    return (protein,aminoseq)
###-------------------------------------------------------------------------------------------------##

#Read csv to dataframe for positive pairs
data_path = os.getcwd() + '/../../DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
colnames = ['DrugBank ID', 'Gene']
df = pd.read_csv(data_path, names = colnames, header = None, index_col = False, sep='\t', skiprows  = [0, 13542])   #rows 0 and 13542 are headers
df['Label'] = 1

#Get the drugs and the smiles
fname = 'drugs_smiles.txt'
path_smiles = os.getcwd() + '/../Data/BIOSNAP/' + fname 
if not os.path.exists(path_smiles):
    get_drug_smiles_drugbank('BIOSNAP')
fname = 'drugs_smiles.txt'
path_smiles = os.getcwd() + '/../Data/BIOSNAP/' + fname 
drugs_smiles = pd.read_csv(path_smiles, delimiter=" ", header = None)
drugs_smiles = dict(zip(drugs_smiles[0], drugs_smiles[1]))
drugs = list(drugs_smiles.keys())

#Obtain the unique genes
genes = df['Gene'].unique()

list_genes=[]
list_targetsequences = []
fname = 'genes_targetsequences.txt'
path_targetsequences = os.getcwd() + '/../Data/BIOSNAP/' + fname 
if not os.path.exists(path_targetsequences):
    for gene in genes:
        gene, targetsequence = getamino_uniprot(gene)
        if gene and targetsequence:
            list_genes.append(gene)
            list_targetsequences.append(targetsequence)
    assert len(list_genes) == len(list_targetsequences),'The length of the Genes  does not match the number of target sequences'
    df_genes = pd.DataFrame()
    df_genes['Gene'] = list_genes
    df_genes['Target Sequence'] = list_targetsequences
    output_path = os.getcwd() + '/../Data/BIOSNAP/genes_targetsequences.txt'
    df_genes.to_csv(output_path, header=None, index = None, sep = ' ')

genes_targetsequences = pd.read_csv(path_targetsequences, delimiter=" ", header = None)
genes_targetsequences = dict(zip(genes_targetsequences[0], genes_targetsequences[1]))
genes = list(genes_targetsequences.keys())  

#Add columns for corresponding SMILES and Target Sequence of each pair
df['SMILES'] = df['DrugBank ID'].map(drugs_smiles)
df['Target Sequence'] = df['Gene'].map(genes_targetsequences)
df = df.dropna()

#Get data for splits
columns = df[['DrugBank ID','Gene', 'Label']]
df_splits = columns.copy()
output_path=os.getcwd() + '/../Data/BIOSNAP/BIOSNAP_pairs.txt'
df_splits.to_csv(output_path)

#Randomnly choose unseen pairs
seenpairs = set(df['DrugBank ID'] + df['Gene'])
unseenpairs = set()
n = len(seenpairs)

df_unseenpairs = pd.DataFrame(columns = ['DrugBank ID', 'Gene', 'Label'])
while n > 0:
    drug_sample = random.choice(drugs)
    gene_sample = random.choice(genes)

    result = drug_sample + gene_sample
    if (result not in seenpairs and result not in unseenpairs):
        sample = pd.DataFrame(data = {'DrugBank ID':[drug_sample], 'Gene': [gene_sample], 'Label':[0]})
        df_unseenpairs = pd.concat([df_unseenpairs, sample], ignore_index = True)
        unseenpairs.add(result)
        n-=1
    
#Add columns for corresponding SMILES and Target Sequence of each unseen pairs
df_unseenpairs['SMILES'] = df_unseenpairs['DrugBank ID'].map(drugs_smiles)
df_unseenpairs['Target Sequence'] = df_unseenpairs['Gene'].map(genes_targetsequences)   

df = pd.concat([df, df_unseenpairs], ignore_index = True)

#Save it as a csv file
output_path = os.getcwd() + '/../Data/BIOSNAP/BIOSNAP.csv'
df.to_csv(output_path)