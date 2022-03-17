import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
from pubchempy import Compound
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
from itertools import repeat

######## SIM MATRIX #########
# SIM MATRIX Calculated via get_drug_sim.py
# read matrix as: (with header & index)
# pd.read_csv(os.path.join(output_path, 'Similarity_Matrix_Drugs.tsv'), sep='\t') 

# RDKIT (default fingerprint, Tanimoto default arguments)
# install: pip3 install rdkit-pypi


'''
Fingerprinting and Molecular Similarity
https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

More details about the algorithm used for the RDKit fingerprint can be found in the “RDKit Book”.

The default set of parameters used by the fingerprinter is: 
	- minimum path size: 1 bond 
	- maximum path size: 7 bonds 
	- fingerprint size: 2048 bits 
	- number of bits set per hash: 2 
	- minimum fingerprint size: 64 bits - target on-bit density 0.0

The default similarity metric used by rdkit.DataStructs.FingerprintSimilarity() 
is the Tanimoto similarity. One can use different similarity metrics:
'''

###########
def get_compound_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles


def get_pairwise_tanimoto(smiles1,smiles2):
	try:
		mol1, mol2 = Chem.MolFromSmiles(str(smiles1)), Chem.MolFromSmiles(str(smiles2))
		fp1, fp2  = Chem.RDKFingerprint(mol1),  Chem.RDKFingerprint(mol2)
		tani = DataStructs.FingerprintSimilarity(fp1,fp2) #pairwise similarity
		return tani
	except:
		return None

def get_smiles_drugbank(drug_entry):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	smiles = None
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
				break
	if not smiles:
		for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
			for ids in exids:
				if(ids.text == 'PubChem Compound'): 
					pubchem_id = exids[1].text
					smiles = get_compound_pubchem(pubchem_id)
					break
	if smiles:
		drug = (drugbank_ID, smiles)
		db_drug_smiles.append(drug)


#############
#############
# change output to (relative path): 
output_path = '../Data/DrugBank'
# read xml
tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
root = tree.getroot()


########### Find Smiles in DrugBank
db_drug_smiles  = []
for drug_entry in tqdm(root):
	get_smiles_drugbank(drug_entry)

# uncomment to save this output into a tsv file in outputfolder
# save file # put in get_drug_sim.py

with open(os.path.join(output_path, 'drugs_w_SMILES.tsv'), 'w') as f:
	_ = f.write('# DrugBank ID\tSMILES\n')
	for item in db_drug_smiles:
		_ = f.write("%s\t%s\n" % item)

len(db_drug_smiles) - len(db_drug_smiles_ok)

db_drug_smiles_ok = []
for sml in tqdm(db_drug_smiles):
	print(sml)
	fp = None
	fp = Chem.MolFromSmiles(str(sml[1]))
	if fp:
		drug = sml
		db_drug_smiles_ok.append(drug)

with open(os.path.join(output_path, 'drugs_w_ok_SMILES.tsv'), 'w') as f:
	_ = f.write('# DrugBank ID\tSMILES\n')
	for item in db_drug_smiles_ok:
		_ = f.write("%s\t%s\n" % item)


####
false_drugs = []
for drug in db_drug_smiles:
	if drug not in db_drug_smiles_ok:
		false_drugs.append(drug[0])

db_drug_smiles_ok[0]
###
# df 
import time

df_smiles = pd.DataFrame(db_drug_smiles_ok, columns=['Drug_ID', 'SMILES'])
df_smiles = df_smiles.set_index('Drug_ID')
df_smiles['sims'] = None  # create an empty column for tanimoto similitude list
## Tanimoto
indxs = df_smiles.index.tolist()  # list of index # can use this later?
all_smiles_ = df_smiles.SMILES.tolist()

all_Sim_Tani = []
for i in range(len(indxs)):
	a = time.time()
	print('='*80)
	print(i+1)
	id_drug_1 = indxs[i]
	smiles_drug_1 = df_smiles.loc[id_drug_1]['SMILES'] 
	tmp = []
	tmp.extend(repeat(smiles_drug_1, len(indxs)) )
	with mp.Pool(processes = mp.cpu_count()-5) as pool:
		results = pool.starmap(get_pairwise_tanimoto, zip(tmp, all_smiles_))
	all_Sim_Tani.append(results)
	b=time.time()
	print(b-a)


df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= indxs, index = indxs) # add columns, & index 

df_all_Sim_Tani.to_csv(os.path.join(output_path, 'Similarity_Matrix_Drugs.tsv'), sep='\t') # add here column %& index next time

##
df_tmp = pd.read_csv('/mnt/md0/data/uveleiro/Similarity_Matrix_Drugs.tsv', index_col=0, sep="\t")

df_tmp.set_index(df_tmp.iloc[0].values)

df_tmp.isna().sum()
df_tmp.iloc[0].isna().sum()


false_drugs[0] in df_tmp.columns.tolist() # returns
df_tmp['DB00515']

df_tmp.dropna(axis=1, how='all')
11297
df_tmp.shape # returns: (11297, 11297)
df_tmp = df_tmp.dropna(axis=1, how='all') # [11297 rows x 11281 columns]
df_tmp = df_tmp.dropna(axis=0, how='all') # [11281 rows x 11281 columns]
if df_tmp.isna().sum().sum() != 0:
	print('Problems, weird NaNs!')
if df_tmp.shape[0] != df_tmp.shape[1]:
	print('Not squared matrix! wot?')