
import argparse
import logging
import os, sys
import shutil
from tqdm import tqdm
import pubchempy as pcp
import subprocess  as sp
import xml.etree.ElementTree as ET

def read_and_extract_DrugBank(path):
	dti_file = os.path.join(path, 'DrugBank_DTIs.tsv')
	drugs = []
	try:
		with open(dti_file, 'r') as f:
			_ = next(f)
			drugs = f.readlines()
			drugs = [drug.strip().split('\t')[0] for drug in tqdm(drugs, desc='Reading DrugBank DB')]
			return drugs
	except FileNotFoundError:
		logging.error('File not found: {}'.format(dti_file))
		sys.exit('Please create the file by running the "get_DTI_from_DrugBank.py" script')

def get_smiles(drug_list, db_path):
	"""
	This function returns the smiles of the drugs in the drug_list.
	"""
	drugs_smiles = []
	dti_file = os.path.join(db_path, 'full_database.xml')
	tree = ET.parse(dti_file)
	root = tree.getroot()
	for drug_entry in tqdm(root, desc='Parsing the database'):
		drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
		logging.debug(f'Drugbank ID: {drugbank_ID}')
		if drugbank_ID in drug_list:
			# get the smile from the DB
			prop = drug_entry.find('{http://www.drugbank.ca}calculated-properties')
			if prop:
				for property in prop:
					if property[0].text.lower() == 'smiles':
						smile = property[1].text
						drugs_smiles.append((drugbank_ID, smile))
						logging.debug(f'{drugbank_ID} - {smile}')
						# continue
			# if not present get it from the drugbank web
			else:
				logging.debug(f'Searching for {drugbank_ID} in pubchem')
				# if not present get it from the pubchem
				drugbank_name = drug_entry.find('{http://www.drugbank.ca}name').text
				smiles =[]
				for comp in pcp.get_compounds(drugbank_name, 'name'):
					smiles.append(comp.canonical_smiles)
				smiles = list(set(smiles))
				if len(smiles) >0:
					drugs_smiles.append((drugbank_ID, smiles[0]))
				else:
					drugs_smiles.append((drugbank_ID, None))
					logging.warning(f'{drugbank_name} not found in anywhere')
	return drugs_smiles

def check_and_create_folder(db_name):
	if not os.path.exists(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/Data', db_name)):
		os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/Data', db_name))

def write_smiles(db_name, drugs_smiles):
	with open(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/Data', db_name, 'drugs_smiles.tsv'), 'w') as f:
		for drug_id, smiles in drugs_smiles:
			_ = f.write(f'{drug_id}\t{smiles}\n')
	return

def create_remove_tmp_folder(path):
	if not os.path.exists(path):
		logging.info('Creating tmp folder: {}'.format(path))
		os.makedirs(path)
		return path
	else: 
		return path

def write_smile_per_file(drugs_smiles, PATH):
	for id, smile in tqdm(drugs_smiles, desc='Writing smiles'):
		if smile:
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
	DB_PATH = args.dbPath
	# DB_PATH =  '/home/margaret/data/jfuente/DTI/Data/DrugBank/'
	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0977-x#Sec19'
	logging.info(f'\n\n{paper_cite}\n')
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = get_DB_name(DB_PATH)

	# read the DB and keep the drugs
	drugs = list(set(read_and_extract_DrugBank(DB_PATH)))
	drugs_smiles = get_smiles(drugs, DB_PATH)
	logging.info(f'{len(drugs_smiles)} drugs found')
	
	write_smiles(db_name, drugs_smiles)

	# write the smiles per file to compute
	PATH = create_remove_tmp_folder(os.path.join('/tmp/SIMCOMP_sucedaneo' , db_name))
	write_smile_per_file(drugs_smiles, PATH)
	# compute the scores
	response = sp.check_call(['java', '-jar', '/home/margaret/data/gserranos/SMILESbasedSimilarityKernels/SMILESSimv2.jar', PATH], )
	if response == 0:
		logging.info('Great succes!\n moving the results to the right place')
		results_folder = os.path.join(os.getcwd() , 'simmatrix')
		destination_folder = '/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/Data' + db_name + '/Drug_simmatrix'
		if not os.path.exists(destination_folder):
			os.makedirs(destination_folder)
		_ = shutil.move(results_folder, destination_folder)
		logging.info(f'The right place is: {destination_folder}')






#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################