import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import json
import multiprocessing as mp
import xml.etree.ElementTree as ET
from itertools import repeat
import logging
import argparse
import argparse
from rdkit import RDLogger                     
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler
import helper_functions_dtinet as hf
from rdkit import Chem

'''
logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)
'''
############### change here for DAVIS !!!!
####
# run with no 0s

#sim
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

def get_protein_sim_matrix(db_name, file_path_SW_pickle, file_path_SW_mat, list_of_protein_nodes, dict_protein_sequence):
	list_of_protein_seqs = [dict_protein_sequence[i] for i in list_of_protein_nodes] # 
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
		with mp.Pool(processes=mp.cpu_count()-20) as pool: # change here for less cpu!
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
	rmtree(tmp_path) # comment for  big matrix ! 


######################################## START MAIN #########################################
#############################################################################################

def main():
	'''
	BindingDB
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

    #######
	### log output detals
	logging.info(
		'''
		This script needs:
			- Common coordinate .tsv files (generated by get_coord.py)
			- DrugBank DTI edge list 
			- smiles dict
			- sequence dict
		Returns:
			- mat*.txt
			- Similarity_Matrix_Drugs.txt
			- Similarity_Matrix_Proteins.txt
		Considering only those nodes that have DTIs
		=> mat_drug_protein.txt not any row/column null
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
	drug_dis = pd.read_csv(os.path.join(wdir, 'edgelist_drug_disease.tsv'), sep='\t')
	drug_se = pd.read_csv(os.path.join(wdir, 'edgelist_drug_se.tsv'), sep='\t')
	drug_drug = pd.read_csv(os.path.join(wdir, 'edgelist_drug_drug.tsv'), sep='\t')
	drug_drug.columns = ['D1', 'D2']
	ppi = pd.read_csv(os.path.join(wdir,'edgelist_PPI.tsv'), sep='\t')
	ppi.columns = ['P1', 'P2']
	prot_dis = pd.read_csv(os.path.join(wdir, 'edgelist_protein_disease.tsv'), sep='\t', usecols=['UniprotID', 'DiseaseID'])
	
	# this is the only one that changes in each model 
	dti = pd.read_csv(os.path.join(wdir, f'DTI_{DB_PATH}.tsv'), sep='\t') 
	dti.columns = ['Drug', 'Protein']
	dti = dti.drop_duplicates()
	logging.info(f'DTI original shape {dti.shape}')

	### SECOND: GET FILTERED NODES 
	logging.info('-'*30)
	logging.info('Gettin Data from tsv files...')

	list_of_drug_nodes = dti.Drug.unique().tolist()
	list_of_protein_nodes = dti.Protein.unique().tolist()
	#wdir = '../Data/Davis_et_al'
	# load data !!! 
	file_path_dic_drug_smiles=  os.path.join(wdir, 'dic_smiles.json')
	file_path_dic_protein_seq = os.path.join(wdir, 'dic_protein_seq.json')

	with open(file_path_dic_drug_smiles, 'r') as f:
		dict_drugid_smiles = json.load(f)

	with open(file_path_dic_protein_seq, 'r') as f:
		dict_protein_sequence = json.load(f)

	# we do not need to take any drop_Drugs now:
	# Before continue, we need to clean from DTI those nodes that are not present in out network
	# these are: list_of_drug_nodes, list_of_protein_nodes
	#drop_drugs = list(set(list_of_drug_nodes).symmetric_difference(dti.Drug.unique().tolist()))
	#drop_proteins = list(set(list_of_protein_nodes).symmetric_difference(dti.Protein.unique().tolist()))
	#index_drop_drugs = dti[dti.Drug.isin(drop_drugs)].index
	#dti = dti.drop(index=index_drop_drugs)
	#index_drop_proteins = dti[dti.Protein.isin(drop_proteins)].index
	#dti = dti.drop(index=index_drop_proteins)
	# now DTI dataframe should not convert into a null row/column
	#list_of_drug_nodes = dti.Drug.unique().tolist()
	#list_of_protein_nodes = dti.Protein.unique().tolist()
	###### do the same for disease and side effect nodes
	# disease nodes for proteins & drugs, but the same
	# diseases linked to drugs
	drop_dis_drugs = list(set(list_of_drug_nodes).symmetric_difference(drug_dis.DrugBankID.unique().tolist()))
	index_drop_drug_di = drug_dis[drug_dis.DrugBankID.isin(drop_dis_drugs)].index
	drug_dis = drug_dis.drop(index=index_drop_drug_di)
	unique_diseases_for_drugs = drug_dis.DiseaseID.unique().tolist()
	# diseases linked to proteins
	drop_dis_protis = list(set(list_of_protein_nodes).symmetric_difference(prot_dis.UniprotID.unique().tolist()))
	index_drop_prot_di = prot_dis[prot_dis.UniprotID.isin(drop_dis_protis)].index
	prot_dis = prot_dis.drop(index=index_drop_prot_di)
	unique_diseases_for_proteins = prot_dis.DiseaseID.unique().tolist()
	# union
	list_of_disease_nodes = list(set(unique_diseases_for_drugs).union(set(unique_diseases_for_proteins)))
	# side effect nodes: drug_se
	drop_se_drugs =  list(set(list_of_drug_nodes).symmetric_difference(drug_se.DrugBank_ID.unique().tolist()))
	index_drop_se = drug_se[drug_se.DrugBank_ID.isin(drop_se_drugs)].index
	drug_se = drug_se.drop(index=index_drop_se)
	list_of_se_nodes = drug_se.se.unique().tolist()

	# Overwrite files for matrix w\o header / index
	np.savetxt(os.path.join(wdir, 'drug.txt'), list_of_drug_nodes, fmt='%s') # list
	np.savetxt(os.path.join(wdir, 'protein.txt'), list_of_protein_nodes, fmt='%s')
	np.savetxt(os.path.join(wdir, 'disease.txt'), list_of_disease_nodes, fmt='%s')
	np.savetxt(os.path.join(wdir, 'se.txt'), list_of_se_nodes, fmt='%s')

	########### 
	# log network information
	logging.info(f'This network has {len(list_of_drug_nodes)} drug nodes')
	logging.info(f'This network has {len(list_of_protein_nodes)} protein nodes')
	### DISEASE NODES 
	logging.info(f'This network has {len(list_of_disease_nodes)} disease nodes')
	### SIDE EFFECT NODES
	logging.info(f'This network has {len(list_of_se_nodes)} side-effect nodes')
	## number of total nodes
	logging.info(f'This network has {len(list_of_protein_nodes) +  len(list_of_drug_nodes) + len(list_of_disease_nodes) + len(list_of_se_nodes)} nodes in total')
	
	########################################################################
	#######################################################################
	
	# THIRD: BUILDING MATRIX
	# once we have the list, we have the index and columns for all matrix!!
	logging.info('-'*30)
	logging.info('Getting matrix....')
	'''
	######### Drug-Disease Matrix
	logging.info('   - Drug Disease Matrix')
	matrix_drug_dis_ = pd.get_dummies(drug_dis.set_index('DrugBankID')['DiseaseID']).max(level=0) 
	#logging.info(f'        * # unique diseases in drugs {len(set(matrix_drug_dis_.columns))} from {len(list_of_disease_nodes)} total nodes')
	matrix_drug_dis = pd.DataFrame(matrix_drug_dis_, columns= list_of_disease_nodes, index= list_of_drug_nodes)
	matrix_drug_dis = matrix_drug_dis.fillna(int(0))
	matrix_drug_dis = matrix_drug_dis.astype(int)
	logging.info(f'        * matrix shape {matrix_drug_dis.shape}')
	logging.info(f'        * # drug-disease assoc edges {matrix_drug_dis.sum().sum()}')
	matrix_drug_dis.to_csv(os.path.join(wdir, 'mat_drug_disease.txt'), index=False, header=False, sep=" ") 

	########## Drug Side Effect Matrix
	logging.info('   - Drug Side Effect Matrix')
	matrix_drug_se_ = pd.get_dummies(drug_se.set_index('DrugBank_ID')['se']).max(level=0)
	#logging.info(f'        * # unique side effects {len(set(matrix_drug_se_.columns))}')
	matrix_drug_se = pd.DataFrame(matrix_drug_se_, columns= list_of_se_nodes, index= list_of_drug_nodes)
	matrix_drug_se = matrix_drug_se.fillna(int(0))
	matrix_drug_se = matrix_drug_se.astype(int)
	logging.info(f'        * Matrix shape {matrix_drug_se.shape}')
	logging.info(f'        * # drug-side effect assoc edges {matrix_drug_se.sum().sum()}')
	matrix_drug_se.to_csv(os.path.join(wdir, 'mat_drug_se.txt'), index=False, header=False, sep=" ") 

	## Drug Drug Matrix  
	logging.info('   - Drug Drug Matrix')
	matrix_drug_drug_ = pd.get_dummies(drug_drug.set_index('D1')['D2']).max(level=0)
	#logging.info(f'        * # unique drugs * that interact {len(set(matrix_drug_drug_.columns))}')
	matrix_drug_drug = pd.DataFrame(matrix_drug_drug_, columns= list_of_drug_nodes, index= list_of_drug_nodes)
	matrix_drug_drug = matrix_drug_drug.fillna(int(0))
	matrix_drug_drug = matrix_drug_drug.astype(int)
	logging.info(f'        * matrix shape {matrix_drug_drug.shape}')
	edges_drug_drug = ( matrix_drug_drug.sum().sum() - np.diag(matrix_drug_drug).sum() )/2 +  np.diag(matrix_drug_drug).sum()
	logging.info(f'        * # drug-drug interaction edges {edges_drug_drug}')
	matrix_drug_drug.to_csv(os.path.join(wdir, 'mat_drug_drug.txt'), index=False, header=False, sep=" ") 

	## Protein-Protein Matrix
	logging.info('   - Protein Protein Matrix')
	ppi_t = ppi[['P2','P1']]
	ppi_t.columns = ['P1', 'P2']
	ppi = ppi.append(ppi_t, ignore_index=True)
	matrix_protein_protein_ = pd.get_dummies(ppi.set_index('P1')['P2']).max(level=0)
	matrix_protein_protein_.columns
	len(set(matrix_protein_protein_.columns))
	#logging.info(f'        * # unique proteins that interact (in HPRD) {len(set(matrix_protein_protein_.columns))}; using prot nodes: {len(list_of_protein_nodes)}')
	matrix_protein_protein = pd.DataFrame(matrix_protein_protein_, columns= list_of_protein_nodes, index= list_of_protein_nodes)
	matrix_protein_protein = matrix_protein_protein.fillna(int(0))
	matrix_protein_protein = matrix_protein_protein.astype(int)
	edges_protein_protein = ( matrix_protein_protein.sum().sum() - np.diag(matrix_protein_protein).sum() )/2 +  np.diag(matrix_protein_protein).sum()
	logging.info(f'        * matrix shape {matrix_protein_protein.shape}')
	logging.info(f'        * # protein-protein edges {edges_protein_protein}')
	matrix_protein_protein.to_csv(os.path.join(wdir, 'mat_protein_protein.txt'), index=False, header=False, sep=" ") 

	# Protein Disease Matrix 
	logging.info('   - Protein Disease Matrix')
	matrix_prot_dis_ = pd.get_dummies(prot_dis.set_index('UniprotID')['DiseaseID']).max(level=0)
	#logging.info(f'        * # unique drugs * that interact {len(set(matrix_prot_dis_.columns))}')
	prot_dis.drop_duplicates()
	matrix_prot_dis = pd.DataFrame(matrix_prot_dis_, columns= list_of_disease_nodes, index= list_of_protein_nodes)
	matrix_prot_dis = matrix_prot_dis.fillna(int(0))
	matrix_prot_dis = matrix_prot_dis.astype(int)
	logging.info(f'        * matrix shape {matrix_prot_dis.shape}')
	logging.info(f'        * # protein-disease edges {matrix_prot_dis.sum().sum()}')
	matrix_prot_dis.to_csv(os.path.join(wdir, 'mat_protein_disease.txt'), index=False, header=False, sep=" ") 
	'''
	# DTI  (DRUG - PROTEIN) ----> Changes for each Database
	logging.info('   - Drug Protein Interaction Matrix (DTIs)')
	dti = dti.drop_duplicates()
	dti.to_csv(os.path.join(wdir, f'final_dtis_{DB_PATH}.tsv'), index=False, header=False, sep="\t")
	matrix_drug_protein_ = pd.get_dummies(dti.set_index('Drug')['Protein']).max(level=0)
	#logging.info(f'        * # unique drugs in DTI info {len(set(matrix_drug_protein_.index))}; # unique drugs in DTI info {len(set(matrix_drug_protein_.columns))}')
	matrix_drug_protein = pd.DataFrame(matrix_drug_protein_, columns= list_of_protein_nodes, index= list_of_drug_nodes)
	matrix_drug_protein = matrix_drug_protein.fillna(int(0))
	matrix_drug_protein = matrix_drug_protein.astype(int)
	logging.info(f'        * matrix shape {matrix_drug_protein.shape}')
	logging.info(f'        * # drug-protein edges {matrix_drug_protein.sum().sum()}')
	'''
	matrix_drug_protein.to_csv(os.path.join(wdir, 'mat_drug_protein.txt'), index=False, header=False, sep=" ") 
	matrix_protein_drug = matrix_drug_protein.T
	matrix_protein_drug.to_csv(os.path.join(wdir, 'mat_protein_drug.txt'), index=False, header=False, sep=" ") 
	
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
	
	################## PROTEIN SIMILARITY MATRIX  
	logging.info('-'*30)
	logging.info('Protein Similarity Matrix....')
	file_path_SW_pickle = os.path.join(wdir, 'prot_sim.pkl')
	file_path_SW_mat = os.path.join(wdir, 'Similarity_Matrix_Proteins.txt')
	if (not os.path.exists(file_path_SW_mat)):
		get_protein_sim_matrix(db_name, file_path_SW_pickle, file_path_SW_mat, list_of_protein_nodes, dict_protein_sequence)
	else:
		logging.info('Matrix already in folder')
	'''
	
################################################################


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################


