import argparse
import logging
import os, sys
import shutil
from re import search
import subprocess  as sp
from tqdm import tqdm
from GetDrugsAndGenes import get_SIMMCOMP_score, get_drug_pubchem
import re


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR']
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

def read_and_extract_BINDING(path):
	pattern = re.compile(r'[\d]+$')
	drugs = []
	with open(path, 'r') as f:
		next(f)
		for line in f:
			line = line.split('\t')
			if line[6] == '':
				drug_id = pattern.search(line[8]).group()
				drug_id = 'BDBM'+drug_id
				drugs.append((drug_id, line[2]))
			else:
				drugs.append((line[6], line[2]))
	return drugs

def check_and_create_folder(db_name):
	if not os.path.exists(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name)):
		os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name))

def write_smiles(db_name, drugs_smiles):
	with open(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name, 'drugs_smiles.tsv'), 'w') as f:
		for drug_id, smiles in drugs_smiles:
			_ = f.write(f'{drug_id}\t{smiles}\n')

def create_remove_tmp_folder(path):
	if not os.path.exists(path):
		logging.info('Creating tmp folder: {}'.format(path))
		os.makedirs(path)
		return path
	else: 
		return path

def write_smile_per_file(drugs_smiles):
	for id, smile in tqdm(drugs_smiles, desc='Writing smiles'):
		with open(os.path.join(PATH, f'{id}.smi'), 'w') as f:
			_ = f.write(smile)

######################################## START MAIN #########################################
#############################################################################################
def main():
	# get the parameters from the user
	parser = argparse.ArgumentParser()
	parser.add_argument("dbPath", help="Path to the database interaction lits",
						default = '/home/margaret/data/jfuente/DTI/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv',
						type=str)
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
						help="Verbosity (between 1-4 occurrences with more leading to more "
							"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
							"DEBUG=4")
	args = parser.parse_args()
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

	# sanity check for the DB
	# DB_PATH = args.dbPath
	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0977-x#Sec19'
	logging.info(f'\n\n{paper_cite}\n')
	DB_PATH =  '/home/margaret/data/jfuente/DTI/Data/BindingDB/BindingDB_All_2021m11/BindingDB_Kd_filtered_columns_subset_HomoSapiens_noduplicities.tsv'
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = get_DB_name(DB_PATH)

	# read the DB and keep the drugs
	drugs = read_and_extract_BINDING(DB_PATH)
	drugs_smiles = list(set(drugs))
	logging.info(f'{len(drugs_smiles)} drugs found')
	
	write_smiles(db_name, drugs_smiles)

	# write the smiles per file to compute
	global PATH
	PATH = create_remove_tmp_folder(os.path.join('/tmp/SIMCOMP_sucedaneo' , db_name))
	write_smile_per_file(drugs_smiles)
	# compute the scores
	response = sp.check_call(['java', '-jar', '/home/margaret/data/gserranos/SMILESbasedSimilarityKernels/SMILESSimv2.jar', PATH], )
	if response == 0:
		logging.info('Great succes!\n moving the results to the right place')
		results_folder = os.path.join(os.getcwd() , 'simmatrix')
		destination_folder = '/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/' + db_name + '/Drug_simmatrix'
		if not os.path.exists(destination_folder):
			os.mkdir(destination_folder)
		_ = shutil.move(results_folder, destination_folder)
		logging.info(f'The right place is: {destination_folder}')






#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################