import os, sys, uuid
import numpy as np
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import multiprocessing as mp
import time
import json
from pubchempy import Compound
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
import xml.etree.ElementTree as ET
from itertools import repeat
import requests
import logging
from re import X, search
import argparse
import argparse
from rdkit import RDLogger                     
import subprocess as sp
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler
import helper_functions_dtinet as hf


# modified 
# only nodes in DTIs!! no null rows

def get_protein_nodes(file_path_prot, file_path_seqs, PPI, prot_dis, DTI):
	#global PPI, prot_dis, DTI
	if (not os.path.exists(file_path_prot)) or (not os.path.exists(file_path_seqs)):
		logging.info('Getting Protein Nodes & Sequences ..')
		''' 
		# PREVIOUS IDEA: 
		## need to do an interesection as well with those that  do not have available sequence
		proteins_linked_to_proteins = set(PPI.P1.tolist() + PPI.P2.tolist()) # 9183 
		proteins_linked_to_disease = set(prot_dis.UniprotID.tolist()) # 8985
		proteins_linked_to_drugs = set(DTI.Protein.tolist()) # 4884 ##### <------------------ THIS IS THE ONLY ONE THAT CHANGES FOR EACH DATABASE
		not_isolated_proteins = list(proteins_linked_to_proteins.union(proteins_linked_to_disease, proteins_linked_to_drugs)) # 11941
		#not_isolated_proteins # not any specific order 
		'''
		'''
		#### SECOND 
		# GET PROTEIN NODES
		# LOAD HPRD NODES & TAKE ONLY THOSE THAT HAVE AT LEAST ONE RELATED DISEASE
		## aqui no estoy teniendo en cuenta las proteinas que estan en DTIs!!! 
		# Load list of HPRD proteins # check that they are in HPRD (node origin)
		proteins_linked_to_proteins = set(PPI.P1.tolist() + PPI.P2.tolist()) # 9183 
		# load protein linked to disease
		proteins_linked_to_disease = set(prot_dis.UniprotID.tolist()) # 8985
		not_isolated_proteins = list(proteins_linked_to_proteins.intersection(proteins_linked_to_disease)) # 8985
		'''
		# TRY ONLY WITH PROTEINS IN DTIs
		# Try only with proteins in dti
		proteins_linked_to_drugs = set(DTI.Protein.tolist()) # 4884 ##### 
		not_isolated_proteins = list(proteins_linked_to_drugs)
		## CHECK SEQUENCE
		logging.info('Checking if sequence is available in Uniprot & retrieving it...')
		list_of_protein_nodes = []
		list_of_protein_seqs  = []
		weird_prots = []
		for prot in tqdm(not_isolated_proteins, desc='Retrieving sequences from Uniprot', position=0, leave=True): 
			prot_id, prot_seq = hf.get_amino_uniprot(prot)
			if(prot_seq):
				list_of_protein_nodes.append(prot_id)
				list_of_protein_seqs.append(prot_seq)
			else:  
				weird_prots.append(prot_id)
		logging.debug(f'The following proteins were not found in uniprot:\n {weird_prots}')
		np.savetxt(os.path.join(file_path_prot), list_of_protein_nodes, fmt='%s')
		with open(file_path_seqs, 'w') as f:
			for i in range(len(list_of_protein_nodes)):
				_ = f.write('>'+list_of_protein_nodes[i]+'\n'+list_of_protein_seqs[i]+'\n')
	else:
		list_of_protein_nodes = np.loadtxt(file_path_prot, dtype='str').tolist()
		list_of_protein_seqs = [ seq for _, seq in list(hf.read_fasta(file_path_seqs)) ]
	#
	return list_of_protein_nodes, list_of_protein_seqs

def get_drug_nodes(file_path_drugs, file_path_dic_drug_smiles, drug_drug, drug_dis, drug_se, DTI):
	# need to define al global the drug matrix !!
	#global drug_drug, drug_dis, drug_se, DTI
	# check
	if ((not os.path.exists(file_path_drugs)) or (not os.path.exists(file_path_dic_drug_smiles))):
		logging.info('Getting Drug Nodes...')
		## get a list of tuples (drug, smiles)
		logging.info('Reading DrugBank xml file')
		tree = ET.parse('../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
		root = tree.getroot()
		logging.info('Retrieving a list of drugs with SMILES')
		# list_drugs_w_smiles  = []
		# for drug_entry in tqdm(root):
		#	list_drugs_w_smiles.append(get_list_drug_w_smiles(drug_entry))
		list_drugs  = []
		list_smiles = []
		for drug_entry in tqdm(root):
			drug_id, smiles = hf.get_smiles(drug_entry)
			if drug_id and smiles:
				list_drugs.append(drug_id)
				list_smiles.append(smiles)
		# 
		assert len(list_drugs) == len(list_smiles), 'The length of the Drug IDs does not match the number of SMILES'
		dict_drugid_smiles = dict(zip(list_drugs, list_smiles))

		'''
		#### PREVIOUS IDEA like in proteins
		# FIND NODES (UNION & INTERSECT) 
		drugs_linked_to_drugs = set(drug_drug.D1.tolist() + drug_drug.D2.tolist()) # 4418
		drugs_linked_to_disease = set(drug_dis.DrugBankID.tolist()) # 2754
		drugs_linked_to_sideeffect = set(drug_se.DrugBank_ID.tolist()) # 675
		drugs_linked_to_proteins = set(DTI.Drug.tolist()) # 7627 ###### <--------- THIS IS THE ONLY ONE THAT CHANGES FOR EACH DATABASE
		not_isolated_drugs = list(drugs_linked_to_drugs.union(drugs_linked_to_disease, drugs_linked_to_sideeffect, drugs_linked_to_proteins )) ### 9720
		not_isolated_drugs.sort()
		'''
		'''
		## NOW TAKE ONLY THOSE DRUGS THAT HAVE DISEASES & THOSE THAT HAVE SE & INTERSECT
		drugs_linked_to_disease = set(drug_dis.DrugBankID.tolist()) # 2754
		drugs_linked_to_sideeffect = set(drug_se.DrugBank_ID.tolist()) # 998
		not_isolated_drugs = drugs_linked_to_disease.intersection(drugs_linked_to_sideeffect) # 860
		'''
		### THIRD TIME
		# Only drugs from DTIs
		drugs_linked_to_proteins = set(DTI.Drug.tolist()) # 7627 ###### <--------- THIS IS THE ONLY ONE THAT CHANGES FOR EACH DATABASE
		not_isolated_drugs = list(drugs_linked_to_proteins)
		# Intersec with available SMILES
		list_of_drug_nodes = list(set(not_isolated_drugs).intersection(set(list_drugs))) #  ##############---->>> ** aqui ya estaria!
		list_of_drug_nodes.sort()
		#len(set(not_isolated_drugs).intersection(set(list_drug_w_smiles))) # 8638
		#
		logging.info(f'There are {len(list_of_drug_nodes)} isolated drugs with available smiles from {len(not_isolated_drugs)} total isolated nodes (diff: {len(not_isolated_drugs)-len(list_of_drug_nodes)})')
		# save files
		logging.info(f'Saving files...')
		np.savetxt(os.path.join(file_path_drugs), list_of_drug_nodes, fmt='%s') # list
		with open(file_path_dic_drug_smiles, 'w', encoding='utf-8') as f:
			json.dump(dict_drugid_smiles, f, ensure_ascii=False, indent=4)
	else:
		logging.info('Reading from existing files...')
		list_of_drug_nodes = np.loadtxt(file_path_drugs, dtype='str').tolist()
		# save json 
		with open(file_path_dic_drug_smiles, 'r') as f:
			dict_drugid_smiles = json.load(f)
	return list_of_drug_nodes, dict_drugid_smiles

#sim
def get_drug_similarity_matrix(file_path_sim_pickle, file_path_simdrug_mat, list_of_drug_nodes, dict_drugid_smiles):
	list_drugs = list_of_drug_nodes
	list_smiles = [dict_drugid_smiles[i] for i in list_drugs] # 
	all_Sim_Tani = []
	dic = {}
	for drug in tqdm(range(len(list_drugs)), desc='Retrieving Pairwise Tanimoto for drugs', position=0, leave=True): 
		id_drug_1 = list_drugs[drug]
		smiles_drug_1 = list_smiles[list_drugs.index(id_drug_1)]
		sim_for_drug = []
		tmp = []
		tmp.extend(repeat(smiles_drug_1, len(list_drugs))) 
		for j in range(len(tmp)):
			result = hf.get_pairwise_tanimoto(tmp[j], list_smiles[j], dic)
			sim_for_drug.append(result)
		all_Sim_Tani.append(sim_for_drug)
	#
	df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= list_drugs, index = list_drugs) # add columns, & index
	logging.info(f'        * Drug Similarity Matrix Shape {df_all_Sim_Tani.shape}')
	# save pickle
	df_all_Sim_Tani.to_pickle(file_path_sim_pickle) # add here column %& index next time
	# save csv for model
	df_all_Sim_Tani.to_csv(file_path_simdrug_mat, sep='\t', header=False, index=False) # add here column %& index next time
# SW
def get_protein_sim_matrix(db_name, file_path_SW_pickle, file_path_SW_mat, list_of_protein_nodes, list_of_protein_seqs):
	targets_seqs = list(zip(list_of_protein_nodes, list_of_protein_seqs))
	# get SW scores
	tmp_path  = hf.create_remove_tmp_folder(os.path.join('/tmp/SmithWaterman' , db_name))
	logging.info(f'Creating temporary folder: {tmp_path}')
	hf.write_all_fastas(targets_seqs, tmp_path)
	all_SmithWaterman = []
	n_targets= len(targets_seqs)
	for pair1 in tqdm(targets_seqs):
		tmp = []
		if not pair1[1]:
			logging.info(f'No sequence for {pair1[0]}')
			continue
		tmp.extend(repeat(pair1, n_targets))
		paths = repeat(tmp_path, n_targets)
		with mp.Pool(processes=mp.cpu_count()-5) as pool:
			results = pool.starmap(hf.get_SW_score, zip(tmp, targets_seqs, paths))
		all_SmithWaterman.append(results)
	
	logging.info('Creating matrix & applying MinMaxScaler() in range (0,1)')
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman, columns=list_of_protein_nodes, index=list_of_protein_nodes)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=list_of_protein_nodes, index=list_of_protein_nodes)
	logging.info(f'        * Protein Similarity Matrix Shape {zscore_SmithWaterman_arr.shape}')
	# save pickle
	zscore_SmithWaterman_arr.to_pickle(file_path_SW_pickle) # add here column %& index next time
	# save csv for model
	zscore_SmithWaterman_arr.to_csv(file_path_SW_mat, sep='\t', header=False, index=False) # add here column %& index next time
	rmtree(tmp_path)



def main():
	'''
	ff
	'''
	parser = argparse.ArgumentParser() 
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
					help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	parser.add_argument('-json', help="If selected, outputs a dictionary in json", action="store_true")
	parser.add_argument("dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)

	args = parser.parse_args()
	RDLogger.DisableLog('rdApp.*')   # disablen warnigns in rdkit
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
	logging.info(
		'''
		This script needs:
			- common coordinate .tsv files (generated by get_coord.py)
			- DTI for each model (as tsv coordinate file)
		
		Returns:
			- mat*.txt
			- Similarity_Matrix_Drugs.txt
			- Similarity_Matrix_Proteins.txt
		'''
		)
	# OUTPUT DIRECTORY
	# sanity check
	DB_PATH = args.dbPath
	logging.info(f'Working in output folder for: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	hf.check_and_create_folder(db_name)
	# Create relative output path
	wdir = os.path.join('../Data', db_name)
	# wdir = '../Data/DrugBank'

	##########################################
	### FIST: LOAD DATA
	logging.info('Loadding Data from tsv files....')
	drug_dis = pd.read_csv(os.path.join(wdir, 'coordinates_drug_disease.tsv'), sep='\t')
	drug_se = pd.read_csv(os.path.join(wdir, 'coordinates_drug_se.tsv'), sep='\t')
	drug_drug = pd.read_csv(os.path.join(wdir, 'coordinates_drug_drug.tsv'), sep='\t')
	drug_drug.columns = ['D1', 'D2']
	PPI = pd.read_csv(os.path.join(wdir,'coordinates_PPI.tsv'), sep='\t')
	PPI.columns = ['P1', 'P2']
	prot_dis = pd.read_csv(os.path.join(wdir, 'coordinates_protein_disease.tsv'), sep='\t', usecols=['UniprotID', 'DiseaseID'])
	# this is the only one that changes in each model ####  
	DTI = pd.read_csv(os.path.join(wdir, f'DTI_{DB_PATH}.tsv'), sep='\t') 
	DTI.columns = ['Drug', 'Protein']
	logging.debug(DTI.head(4))
	DTI = DTI.drop_duplicates()

	### SECOND: GET FILTERED NODES 
	logging.info('-'*30)
	logging.info('Gettin Data from tsv files...')
	# DRUG NODES
	# file_path_drugs = os.path.join(wdir, 'tmp_drugs.txt')
	file_path_drugs = os.path.join(wdir, 'drug.txt')
	file_path_dic_drug_smiles=  os.path.join(wdir, 'dic_smiles.json')
	list_of_drug_nodes, dict_drugid_smiles = get_drug_nodes(file_path_drugs, file_path_dic_drug_smiles, drug_drug, drug_dis, drug_se, DTI)
	logging.info(f'This network has {len(list_of_drug_nodes)} drug nodes')
	# PROTEIN NODES
	file_path_prot = os.path.join(wdir, 'protein.txt')
	file_path_seqs = os.path.join(wdir, 'protein_seqs.fasta')
	list_of_protein_nodes, list_of_protein_seqs = get_protein_nodes(file_path_prot, file_path_seqs, PPI, prot_dis, DTI)
	logging.info(f'This network has {len(list_of_protein_nodes)} protein nodes')
	# DISEASE NODES 
	list_of_disease_nodes = set(drug_dis.DiseaseID).union(set(prot_dis.DiseaseID)) ## union?
	logging.info(f'This network has {len(list_of_disease_nodes)} disease nodes')
	# SIDE EFFECT NODES
	list_of_se_nodes = drug_se.se.unique().tolist()
	logging.info(f'This network has {len(list_of_se_nodes)} side-effect nodes')
	## number of total nodes
	logging.info(f'This network has {len(list_of_protein_nodes) +  len(list_of_drug_nodes) + len(list_of_disease_nodes) + len(list_of_se_nodes)} nodes in total')

	'''
	#################
	# THIRD: BUILDING MATRIX
	# once we have the list, we have the index and columns for all matrix!!
	logging.info('-'*30)
	logging.info('Getting matrix....')

	######### Drug-Disease Matrix
	logging.info('   - Drug Disease Matrix')
	matrix_drug_dis_ = pd.get_dummies(drug_dis.set_index('DrugBankID')['DiseaseID']).max(level=0) 
	logging.info(f'        * # unique diseases in drugs {len(set(matrix_drug_dis_.columns))} from {len(list_of_disease_nodes)} total nodes')
	matrix_drug_dis = pd.DataFrame(matrix_drug_dis_, columns= list_of_disease_nodes, index= list_of_drug_nodes)
	matrix_drug_dis = matrix_drug_dis.fillna(int(0))
	matrix_drug_dis = matrix_drug_dis.astype(int)
	logging.info(f'        * matrix shape {matrix_drug_dis.shape}')
	logging.info(f'        * # drug-disease assoc edges {matrix_drug_dis.sum().sum()}')
	matrix_drug_dis.to_csv(os.path.join(wdir, 'mat_drug_disease.txt'), index=False, header=False, sep=" ") 

	########## Drug Side Effect Matrix
	logging.info('   - Drug Side Effect Matrix')
	matrix_drug_se_ = pd.get_dummies(drug_se.set_index('DrugBank_ID')['se']).max(level=0)
	logging.info(f'        * # unique side effects {len(set(matrix_drug_se_.columns))}')
	matrix_drug_se = pd.DataFrame(matrix_drug_se_, columns= matrix_drug_se_.columns, index= list_of_drug_nodes)
	matrix_drug_se = matrix_drug_se.fillna(int(0))
	matrix_drug_se = matrix_drug_se.astype(int)
	logging.info(f'        * Matrix shape {matrix_drug_se.shape}')
	logging.info(f'        * # drug-side effect assoc edges {matrix_drug_se.sum().sum()}')
	matrix_drug_se.to_csv(os.path.join(wdir, 'mat_drug_se.txt'), index=False, header=False, sep=" ") 

	## Drug Drug Matrix  
	logging.info('   - Drug Drug Matrix')
	matrix_drug_drug_ = pd.get_dummies(drug_drug.set_index('D1')['D2']).max(level=0)
	logging.info(f'        * # unique drugs * that interact {len(set(matrix_drug_drug_.columns))}')
	matrix_drug_drug = pd.DataFrame(matrix_drug_drug_, columns= list_of_drug_nodes, index= list_of_drug_nodes)
	matrix_drug_drug = matrix_drug_drug.fillna(int(0))
	matrix_drug_drug = matrix_drug_drug.astype(int)
	logging.info(f'        * matrix shape {matrix_drug_drug.shape}')
	logging.info(f'        * # drug-drug interaction edges {matrix_drug_drug.sum().sum()/2}')
	matrix_drug_drug.to_csv(os.path.join(wdir, 'mat_drug_drug.txt'), index=False, header=False, sep=" ") 

	## Protein-Protein  Matrix
	logging.info('   - Protein Protein Matrix')
	PPI_t = PPI[['P2','P1']]
	PPI_t.columns = ['P1', 'P2']
	PPI = PPI.append(PPI_t, ignore_index=True)
	matrix_protein_protein_ = pd.get_dummies(PPI.set_index('P1')['P2']).max(level=0)
	matrix_protein_protein_.columns
	len(set(matrix_protein_protein_.columns))
	logging.info(f'        * # unique proteins that interact (in HPRD) {len(set(matrix_protein_protein_.columns))}; using prot nodes: {len(list_of_protein_nodes)}')
	matrix_protein_protein = pd.DataFrame(matrix_protein_protein_, columns= list_of_protein_nodes, index= list_of_protein_nodes)
	matrix_protein_protein = matrix_protein_protein.fillna(int(0))
	matrix_protein_protein = matrix_protein_protein.astype(int)
	logging.info(f'        * matrix shape {matrix_protein_protein.shape}')
	logging.info(f'        * # protein-protein edges {matrix_protein_protein.sum().sum()/2}')
	matrix_protein_protein.to_csv(os.path.join(wdir, 'mat_protein_protein.txt'), index=False, header=False, sep=" ") 

	# Protein Disease Matrix 
	logging.info('   - Protein Disease Matrix')
	matrix_prot_dis_ = pd.get_dummies(prot_dis.set_index('UniprotID')['DiseaseID']).max(level=0)
	logging.info(f'        * # unique drugs * that interact {len(set(matrix_prot_dis_.columns))}')
	prot_dis.drop_duplicates()
	matrix_prot_dis = pd.DataFrame(matrix_prot_dis_, columns= list_of_disease_nodes, index= list_of_protein_nodes)
	matrix_prot_dis = matrix_prot_dis.fillna(int(0))
	matrix_prot_dis = matrix_prot_dis.astype(int)
	logging.info(f'        * matrix shape {matrix_prot_dis.shape}')
	logging.info(f'        * # protein-disease edges {matrix_prot_dis.sum().sum()}')
	matrix_prot_dis.to_csv(os.path.join(wdir, 'mat_protein_disease.txt'), index=False, header=False, sep=" ") 

	# DTI  (DRUG - PROTEIN) ----> Changes for each Database
	logging.info('   - Drug Protein Interaction Matrix (DTIs)')
	DTI = DTI.drop_duplicates()
	matrix_drug_protein_ = pd.get_dummies(DTI.set_index('Drug')['Protein']).max(level=0)
	logging.info(f'        * # unique drugs in DTI info {len(set(matrix_drug_protein_.index))}; # unique drugs in DTI info {len(set(matrix_drug_protein_.columns))}')
	matrix_drug_protein = pd.DataFrame(matrix_drug_protein_, columns= list_of_protein_nodes, index= list_of_drug_nodes)
	matrix_drug_protein = matrix_drug_protein.fillna(int(0))
	matrix_drug_protein = matrix_drug_protein.astype(int)
	logging.info(f'        * matrix shape {matrix_drug_protein.shape}')
	logging.info(f'        * # drug-protein edges {matrix_drug_protein.sum().sum()}')
	# pd.Series(np.diag(df), index=[df.index, df.columns]) ejemplo de contar 1s en la diagonal!! cambiar aqui
	matrix_drug_protein.to_csv(os.path.join(wdir, 'mat_drug_protein.txt'), index=False, header=False, sep=" ") 
	'''

	# Drug Similarity Matrix
	logging.info('-'*30)
	logging.info('Drug Similarity matrix....')
	# definir file path
	# call function directly
	file_path_sim_pickle = os.path.join(wdir, 'drug_sim.pkl')
	file_path_simdrug_mat = os.path.join(wdir, 'Similarity_Matrix_Drugs.txt')
	if ((not os.path.exists(file_path_simdrug_mat)) ):
		# call function: 
		logging.info('Calculating matrix...')
		get_drug_similarity_matrix(file_path_sim_pickle, file_path_simdrug_mat, list_of_drug_nodes, dict_drugid_smiles)
		#pass
	else:
		logging.info('Matrix already in folder')
	'''
	################## PROTEIN SIMILARITY MATRIX  
	logging.info('-'*30)
	logging.info('Protein Similarity Matrix....')
	file_path_SW_pickle = os.path.join(wdir, 'prot_sim.pkl')
	file_path_SW_mat = os.path.join(wdir, 'Similarity_Matrix_Proteins.txt')
	if (not os.path.exists(file_path_SW_mat)):
		get_protein_sim_matrix(db_name, file_path_SW_pickle, file_path_SW_mat, list_of_protein_nodes, list_of_protein_seqs)
	else:
		logging.info('Matrix already in folder')
	'''	
		

#############################################



#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
