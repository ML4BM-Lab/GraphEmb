import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
from pubchempy import Compound
import json
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
		return np.nan

#############

# change output to (relative path): 
output_path = '../Data/DrugBank'
# read xml
tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
root = tree.getroot()


########### Find Smiles in DrugBank
db_drug_smiles  = []
for drug_entry in tqdm(root):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	smiles = np.nan
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
				break
	if smiles == np.nan:
		for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
			for ids in exids:
				if(ids.text == 'PubChem Compound'): 
					pubchem_id = exids[1].text
					smiles = get_compound_pubchem(pubchem_id)
					break
	drug = (drugbank_ID, smiles)
	db_drug_smiles.append(drug)

### uncomment to save this output into a tsv file in outputfolder
# save file # put in get_drug_sim.py
# with open(os.path.join(output_path, 'drugs_w_SMILES.tsv'), 'w') as f:
#	#_ = f.write('# DrugBank ID\tProtein list\n')
#	for item in db_drug_smiles:
#		_ = f.write("%s\t%s\n" % item)
 
t = pd.DataFrame(db_drug_smiles, columns=['Drug', 'SMILES']).dropna()
list_drug_w_smiles = t.Drug.tolist()
np.savetxt(os.path.join(output_path ,'drug_w_smiles.txt'),  list_drug_w_smiles  , newline='\n', fmt='%s')


# df 
df_smiles = pd.DataFrame(db_drug_smiles, columns=['Drug_ID', 'SMILES'])
df_smiles = df_smiles.set_index('Drug_ID')
df_smiles = df_smiles.dropna()
df_smiles['sims'] = np.nan  # create an empty column for tanimoto similitude list

## Tanimoto
indxs = df_smiles.index.tolist()  # list of index # can use this later?
all_smiles_ = df_smiles.SMILES.tolist()
all_Sim_Tani = []
for i in range(len( indxs )):
	print('='*80)
	print(i+1)
	id_drug_1 = indxs[i]
	smiles_drug_1 = df_smiles.loc[id_drug_1]['SMILES'] 
	tmp = []
	tmp.extend(repeat(smiles_drug_1, len(indxs)) )
	with mp.Pool(processes = mp.cpu_count()-5) as pool:
		results = pool.starmap(get_pairwise_tanimoto, zip(tmp, all_smiles_))
	all_Sim_Tani.append(results)


df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= indxs, index = indxs) # add columns, & index

df_all_Sim_Tani.to_csv(os.path.join(output_path, 'Similarity_Matrix_Drugs.tsv'), sep='\t')
