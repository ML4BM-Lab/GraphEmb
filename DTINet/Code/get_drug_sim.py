import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd
import numpy as np
import os,sys
from pubchempy import Compound
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
from itertools import repeat
import argparse
import logging
from re import search


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


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank', 'E', 'GPCR', 'IC', 'NR']
	for db in DB_NAMES:
		if search(db, path):
			logging.info(f'Database: {db}')
			if db in ['E', 'GPCR', 'IC', 'NR']:
				db = os.path.join('Yamanashi_et_al_GoldStandard', db)
				return db
			else:
				return db
	logging.error(f'Database: {db} not found')
	sys.exit('Please provide a valid database')


def check_and_create_folder(db_name):
	if not os.path.exists(os.path.join('../Data', db_name)):
		os.mkdir(os.path.join('../Data', db_name))


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


def get_smiles_drugbank(drug_entry, db_drug_smiles):
	'''
	Get smiles from DrugBank xml, 
	for SMILES that are not in DrugBank retrieves them from PubChem
	Then, check if we can create a fingerprint (if not, do not include)
	# Note: without fp check: len == 11310
	# with fp check: len == 11294  # --> 16 disappear! 
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
	if smiles:
		# check if fp exists
		fp = Chem.MolFromSmiles(str(smiles))
		if fp:
			drug = (drugbank_ID, smiles)
			db_drug_smiles.append(drug)





######################################## START MAIN #########################################
#############################################################################################

def main():
	'''
	CTD
	'''
	parser = argparse.ArgumentParser() 
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
					help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	parser.add_argument("-o", help="Output path ", default='../Data', type=str)
	args = parser.parse_args()
	# -) log info; 
	# define logging level
	log_levels = {
		0: logging.CRITICAL,
		1: logging.ERROR,
		2: logging.WARN,
		3: logging.INFO,
		4: logging.DEBUG,
	}
	# set the logging info
	level= log_levels[args.verbosity]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)

	#######
    ###  output details
	logging.info(
		'''
		This script generates a Similarity Matrix of Drugs 
		using SMILES and Fingerprints + Tanimoto (with RDKit)
		For DTINet !
		'''
		)
	##### PATHS ???? ##### ---------------------------> *****
	output_path = args.oPath

	## CTD DATA FOLDER
	# info 
	logging.info(f'Processing SMILES in {output_path}')

	logging.info('Reading xml file...')
	tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/DrugBank/full_database.xml')
	root = tree.getroot()

	## get a list of tuples (drug, smiles)
	db_drug_smiles  = []
	for drug_entry in tqdm(root):
		get_smiles_drugbank(drug_entry, db_drug_smiles)
	print(len(db_drug_smiles))

	# process data
	df_smiles = pd.DataFrame(db_drug_smiles, columns=['Drug_ID', 'SMILES'])
	df_smiles = df_smiles.set_index('Drug_ID')
	df_smiles['sims'] = None  # create an empty column for tanimoto similitude list
	## Tanimoto
	indxs = df_smiles.index.tolist()  # list of index # can use this later?
	all_smiles_ = df_smiles.SMILES.tolist()

	all_Sim_Tani = []
	all_Sim_Tani = []
	for i in range(len(indxs[:5])): ### CHANGE LATER
		print('='*80)
		print(i+1)
		id_drug_1 = indxs[i]
		smiles_drug_1 = df_smiles.loc[id_drug_1]['SMILES'] 
		tmp = []
		tmp.extend(repeat(smiles_drug_1, len(indxs[:5]))) ### CHANGE LATER!!!
		with mp.Pool(processes = mp.cpu_count()-5) as pool:
			results = pool.starmap(get_pairwise_tanimoto, zip(tmp, all_smiles_))
		all_Sim_Tani.append(results)
	## CHANGE LATER
	df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= indxs[:5], index = indxs[:5]) # add columns, & index 
	df_all_Sim_Tani.to_csv(os.path.join(output_path, 'test_matrix.tsv'), sep='\t', header=True, index=True) # add here column %& index next time


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################

