import argparse
import logging
import os
import shutil
import subprocess  as sp
import re
import helper_functions_DTI2Vec as hf


def BINDINGDB_SIMCOMP(DB_PATH =  './../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'):

	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	JAR_PATH = '/home/margaret/data/gserranos/SMILESbasedSimilarityKernels/SMILESSimv2.jar'

	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0977-x#Sec19'
	logging.info(f'\n\n{paper_cite}\n')
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	# read the DB and keep the drugs
	drugs_smiles = hf.read_and_extract_BINDING_smiles(DB_PATH)
	logging.info(f'{len(drugs_smiles)} drugs found')
	
	# write the smiles per file to compute
	tmp_path = hf.write_smile_per_file(drugs_smiles, db_name)
	# compute the scores
	function_name = 'TFIDF'
	logging.info(f'Computing the scores\n This can take several minutes...')
	response = sp.check_call(['java', '-jar', JAR_PATH, tmp_path, function_name])
	if response == 0:
		logging.info('Great succes!\n moving the results to the right place')
		results_folder = os.path.join(os.getcwd() , 'simmatrix')
		destination_folder = './../Data/' + db_name + '/Drug_simmatrix'
		if not os.path.exists(destination_folder):
			os.mkdir(destination_folder)
		_ = shutil.move(results_folder, destination_folder)
		logging.info(f'The right place is: {destination_folder}')

def YAMANASHI_SIMCOMP():
	pass

def DAVIS_SIMCOMP():
	pass

def BIOSNAP_SIMCOMP():
	pass

def DRUGBANK_SIMCOMP():
	pass

