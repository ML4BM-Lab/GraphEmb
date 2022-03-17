import numpy as np
import pandas as pd
import random
import requests
import os
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
    return aminoseq
###-------------------------------------------------------------------------------------------------##

#Read csv to dataframe for positive pairs
data_path = os.getcwd() + '/../../../Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'

#df = pd.read_csv('ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv', index_col = False, sep='\t', skiprows  = [13542])   #row 13542 is a header
colnames = ['DrugBank ID', 'Gene']
df = pd.read_csv(data_path, names = colnames, header = None, index_col = False, sep='\t', skiprows  = [0, 13542])   #rows 0 and 13542 are headers
df['Label'] = [1]*len(df.index)
#print(df)

#Obtain the unique drugs and genes
drugs = df['DrugBank ID'].unique()
genes = df['Gene'].unique()

#print(drugs)
#print(genes)

#Get the smiles and sequences of drugs and proteins

fname = 'smiles_BIOSNAP.txt'
path_smiles = os.getcwd() + '/' + fname
    
if not os.path.exists(path_smiles):
    smiles = np.vectorize(getdrug_drugbank)(drugs)
    np.savetxt(fname, smiles, fmt="%s")
else:
    smiles = np.genfromtxt(fname, dtype = 'str')

fname = 'targetsequences_BIOSNAP.txt'
path_targetsequences = os.getcwd() + '/' + fname
    
if not os.path.exists(path_targetsequences):
    targetsequences = np.vectorize(getamino_uniprot)(genes)
    np.savetxt(fname, targetsequences, fmt="%s")
    
else:
    targetsequences = np.genfromtxt(fname, dtype = 'str')

#print(smiles)
#print(targetsequences)

#Create dictionaries
drugs_smiles = {d:s for d, s in zip(drugs, smiles)}
genes_targetsequences = {g:t for g,t in zip(genes, targetsequences)}

#Randomnly choose unseen pairs
seenpairs = set(df['DrugBank ID'] + df['Gene']) #create set of strings of drug+gene
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
output_path = os.getcwd() + '/../Data/BIOSNAP/BIOSNAP.csv'
df.to_csv(output_path)