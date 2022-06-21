import os
import numpy as np
import pandas as pd
import logging 
import glob
from soupsieve import select
from tqdm import tqdm
import os

from rdkit import Chem


def get_interactions(DTIs):
    Proteins = DTIs.Protein.unique().tolist()
    Drugs = DTIs.Drug.unique().tolist()
    interactions = []
    for drug in tqdm(Drugs):
        drug_list = []
        drug_list.append(drug)
        # positive links
        positives =  DTIs.Protein[DTIs.Drug == drug].tolist()
        drug_list.append(positives)
        negatives = list(set(Proteins).difference(set(positives)))
        drug_list.append(negatives)
        # finish list
        interactions.append(drug_list)
    return interactions


DTIs = pd.read_csv('../../DTINet/Data/Yamanashi_et_al_GoldStandard/NR/final_dtis_NR.tsv', sep = "\t", names=['Drug', 'Protein'])
#DTIs = pd.read_csv('../../DB/Data/DrugBank/DrugBank_DTIs.tsv', sep = "\t")
#DTIs.columns = ['Drug', 'Protein']
#DTIs = pd.read_csv('../../DTINet/Data/Yamanashi_et_al_GoldStandard/E/final_dtis_E.tsv', sep = "\t", names=['Drug', 'Protein'])


prots = DTIs.Protein.unique().tolist()
#######
prot_pdb = glob.glob('../Data/Clean_from_PDB/*.pdb')
prot_pdb = [prot.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for prot in prot_pdb]
prot_alpha = glob.glob('../Data/Clean_from_AFold/*.pdb')
prot_alpha = [prot.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for prot in prot_alpha]
aval_prots = prot_pdb + prot_alpha
prots_w_str = set(prots).intersection(set(aval_prots))
logging.info(f'Avaliable proteins with structure: {len(prots_w_str)}')
# 
e_in_pdb = list(prots_w_str.intersection(set(prot_pdb)))
e_in_afold = list(prots_w_str.intersection(set(prot_alpha)))


#np.savetxt('../Data/test_E_docking/test_E_docking_frompdb.txt', e_in_pdb, fmt='%s')
#np.savetxt('../Data/test_E_docking/test_E_docking_fromafold.txt', e_in_afold, fmt='%s')


np.savetxt('../Data/test_NR_docking/test_NR_docking_frompdb.txt', e_in_pdb, fmt='%s')
np.savetxt('../Data/test_NR_docking/test_NR_docking_fromafold.txt', e_in_afold, fmt='%s')

# yamanishi same type of proteins
# let us take a lower threshold zB 2
# get matrix or from matrix
#
# load rmsd data
data_rmsd = pd.read_csv('../Data/matrix_rmsd_NR.tsv', sep="\t")
list_prot = [colname.replace('_A', '') for colname in data_rmsd.columns.tolist()]
data_rmsd.columns, data_rmsd.index = list_prot, list_prot
data_rmsd

# thresholds
threshold1 =  data_rmsd.mean().mean() #2.5
mask_data = data_rmsd < threshold1

# DTIs with available structure for targets
DTIs_structural = DTIs[DTIs.Protein.isin(prots_w_str)]
interactions = get_interactions(DTIs_structural)

# statistics

uniques = (len(DTIs_structural.Drug.unique()), len(DTIs_structural.Protein.unique()))

n_pos = len(DTIs_structural)
n_neg = uniques[0]*uniques[1]
n_pos
#
plausible_negatives = []
for i in interactions:
    drug_id = i[0]
    positive = i[1]
    negatives = i[2]
    if len(positive) == 1: # do a loop here and take all if docking? instead of this? and [d1][][], [d1][][]
        result_ = data_rmsd[mask_data][positive].loc[negatives].dropna().sort_values(by=positive).reset_index()
        result_.columns = ['uniprotid', 'rmsd']
        # list of tuples of plausible negatives
        plausible_negatives_4drug = [(drug_id, uni) for uni in result_.uniprotid] 
        plausible_negatives.extend(plausible_negatives_4drug)
    else:
        print(f'for {drug_id} more than one positive protein')


print(f'Current balance. Positives: {n_pos}, Negatives: {len(plausible_negatives)}')


##########################################
##########################################
### use this loop
### try all with nr
plausible_negatives = []
for i in interactions:
    drug_id, positives, negatives = i[0], i[1], i[2]
    #print(drug_id)
    for positive in positives:
        #if len(positive) == 1: # do a loop here and take all if docking? instead of this? and [d1][][], [d1][][]
        #positive
        result_ = data_rmsd[mask_data][positive].loc[negatives].dropna().reset_index().sort_values(by=positive)
        result_.columns = ['uniprotid', 'rmsd']
        # list of tuples of plausible negatives
        plausible_negatives_4pos = [(drug_id, uni) for uni in result_.uniprotid] 
        plausible_negatives.extend(plausible_negatives_4pos)
        #if len(positives) > 1:
        #    print(f'for {drug_id} more than one positive protein')
        #break



#### random selection of 10 positive DTI PAIRS and 10 negative pairs
import random

true_positives = DTIs_structural.to_records(index=False).tolist()

random.seed(5)
selected_positives = random.sample(true_positives, 10)
selected_negatives = random.sample(plausible_negatives, 10)

np.savetxt('../Data/moldock_10pairs/pairs_pos.txt', selected_positives, fmt='%s')
np.savetxt('../Data/moldock_10pairs/pairs_neg.txt', selected_negatives, fmt='%s')

# 
selected_positives = np.loadtxt('../Data/moldock_10pairs/pairs_pos.txt', dtype=str)
selected_negatives = np.loadtxt('../Data/moldock_10pairs/pairs_neg.txt', dtype=str)

len(set(selected_positives[:,1].tolist()).union(set(selected_negatives[:,1].tolist())))

drugs = list(set(selected_positives[:,0].tolist()).union(set(selected_negatives[:,0].tolist())))
drugs


import xml.etree.ElementTree as ET
from pubchempy import Compound


def get_compound_pubchem(drug):
	return Compound.from_cid(drug).isomeric_smiles

def get_smiles(drug_entry):
	'''
	Get list of drugs with smiles from DrugBank xml, 
	for SMILES that are not in DrugBank retrieves them from PubChem
	Then, check if we can create a fingerprint (if not, do not include)
	'''
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	smiles = None
	fp = None
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
	if not smiles:
		return(drugbank_ID, None)
	elif Chem.MolFromSmiles(str(smiles)):
		return(drugbank_ID, smiles)
	return(drugbank_ID, None)



tree = ET.parse(os.path.join(os.getcwd(), '../../DB/Data/DrugBank/full_database.xml'))
root = tree.getroot()
logging.info('Retrieving a list of drugs with SMILES')

list_drugs  = []
list_smiles = []
for drug_entry in tqdm(root, position=0, leave=True):
    drug_id, smiles = get_smiles(drug_entry)
    if drug_id and smiles:
        list_drugs.append(drug_id)
        list_smiles.append(smiles)
# 
assert len(list_drugs) == len(list_smiles), 'The length of the Drug IDs does not match the number of SMILES'
dict_drugid_smiles = dict(zip(list_drugs, list_smiles))




for drug in drugs:
    mol = Chem.MolFromSmiles(dict_drugid_smiles[drug])
    mol_H = Chem.AddHs(mol)
    #print(Chem.MolToMolBlock(mol_H))
    path = f'../Data/moldock_10pairs/ligands_sdf/{drug}.sdf'
    with Chem.SDWriter(path) as w:
        w.write(mol_H)




fix_protein(filename='1AZ8_clean.pdb',addHs_pH=7.4,try_renumberResidues=True,output='1AZ8_clean_H.pdb')


# need permission to mkdir in /share! 
# https://apolo-docs.readthedocs.io/en/latest/software/applications/lePro/2013/index.html

