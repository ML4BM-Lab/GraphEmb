import numpy as np
import pandas as pd
import random
import requests
import os
import xml.etree.ElementTree as ET
random.seed(1)

###----------------------------------DRUGBANK FDA-------------------------------------###
#The drugs can be obtained from the link 
#https://go.drugbank.com/structures/small_molecule_drugs/{DRUG}.smiles
#Example -> Drug: DB00755 -> https://go.drugbank.com/structures/small_molecule_drugs/DB00755.smiles

#The proteins can be obtained from the uniprot link
#https://www.uniprot.org/uniprot/{PROTEIN}.fasta
#Example -> Protein: Q16539 -> https://www.uniprot.org/uniprot/Q16539.fasta

def get_compound_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles


def getdrug_drugbank(drug):
    r = requests.get(f'https://go.drugbank.com/structures/small_molecule_drugs/{drug}.smiles')
    if r.text=='':
        for exids in drug.findall('.//{http://www.drugbank.ca}external-identifier'):
            for ids in exids:
                if (ids.text == 'PubChem Compound'):
                    pubchem_id = exids[1].text
                    smiles = get_compound_pubchem(pubchem_id)
                    print(smiles)
                    return smiles
    return r.text

	

def getamino_uniprot(protein):
    r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

###----------------------------------DRUGBANK FDA-------------------------------------###

#Read csv to dataframe for positive pairs
data_path = os.getcwd() + '/../../../Data/DrugBank/DrugBank_DTIs.tsv'
colnames = ['DrugBank ID', 'Gene']
df = pd.read_csv(data_path, names = colnames, index_col=False, sep='\t', skiprows = [0])
df['Label'] = [1]*len(df.index)
print(df)

#Obtain the unique drugs and genes
drugs = df['DrugBank ID'].unique()
genes = df['Gene'].unique()

#print(drugs)
#print(genes)

#Get the smiles and sequences of drugs and proteins
fname = 'smiles_DrugBank.txt'
path_smiles = os.getcwd() + '/../Data/DrugBank/' + fname 
if not os.path.exists(path_smiles):
    tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/DrugBank/full_database.xml')
    root = tree.getroot()
    smiles  = []
    for drug_entry in tqdm(root):
        print(drug_entry)
		#smiles.append(get_list_drug_w_smiles(drug_entry))
    smiles = np.vectorize(getdrug_drugbank)(drugs)
    np.savetxt(path_smiles, smiles, fmt="%s")
else:
    smiles = np.genfromtxt(path_smiles, dtype = 'str')

fname = 'targetsequences_DrugBank.txt'
path_targetsequences = os.getcwd() + '/../Data/DrugBank/' + fname
    
if not os.path.exists(path_targetsequences):
    targetsequences = np.vectorize(getamino_uniprot)(genes)
    np.savetxt(path_targetsequences, targetsequences, fmt="%s")
    
else:
    targetsequences = np.genfromtxt(path_targetsequences, dtype = 'str')

#print(smiles)
#print(targetsequences)

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

cols = ['DrugBank ID', 'Gene', 'SMILES', 'Target Sequence', 'Label']
df = df[cols]

#Save it as a csv file
output_path = os.getcwd() + '/../Data/DrugBank/DrugBank.txt'
df.to_csv(output_path, header=None, index = None, sep = ' ', na_rep = ' ')