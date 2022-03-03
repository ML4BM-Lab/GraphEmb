import argparse
import logging
import os, sys, uuid
import requests
import pandas as pd
from re import search
from tqdm import tqdm
from itertools import repeat
import subprocess as sp
import multiprocessing as mp
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler


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
				try:
					os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db))
					return db
				except FileExistsError:
					return db
			else:
				try:
					os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db))
					return db
				except FileExistsError:
					return db
	logging.error(f'Database: {db} not found')
	sys.exit('Please provide a valid database')

def read_and_extract_targets(path):
	"""
	This function reads the database and returns the targets.
	"""
	targets = []
	with open(path, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				_, target = line.split('\t')
				targets.append(target.strip())
	return targets


def extract_score(results_file):
	with open(results_file, 'r') as f:
		for line in f:
			if not line.startswith('# Score:'):
				continue
			else:
				return float(line.split()[-1])

def create_remove_tmp_folder(path):
	if not os.path.exists(path):
		logging.info('Creating tmp folder: {}'.format(path))
		os.makedirs(path)
		return path
	else: 
		return path

def write_fasta(path, target, seq):
	fl_name = os.path.join(path, target.replace(':', '_')+'.fasta')
	if os.path.exists(fl_name):
		logging.debug(f'File {fl_name} already exists')
		return fl_name
	with open(fl_name, 'w') as f:
		_ = f.write('>'+target+'\n'+seq+'\n')
	return fl_name

def check_and_create_fasta(target, seq):
	global PATH
	fasta1 = os.path.join(PATH, target.replace(':', '_')+'.fasta')
	if not os.path.exists(fasta1):
		fasta1 = write_fasta(PATH, target, seq)
	return fasta1

def get_SW_score(pair1, pair2):
	global PATH
	target1, seq1 = pair1
	target2, seq2 = pair2
	fasta1 = os.path.join(PATH, target1.replace(':', '_')+'.fasta')
	if not os.path.exists(fasta1):
		fasta1 = write_fasta(PATH, target1, seq1)
	fasta2 = os.path.join(PATH, target2.replace(':', '_')+'.fasta')
	fasta2 = write_fasta(PATH, target2, seq2)
	result_ID = str(uuid.uuid4())
	result_file = os.path.join(PATH, result_ID+'_results.txt')
	args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
			'-asequence', fasta1 , '-bsequence', fasta2, 
			'-gapopen', '10.0', '-gapext', '0.5', 
			'-outfile', result_file]
	try:
		_ = sp.check_call(args, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		score = extract_score(result_file)
		os.remove(result_file)
		return score
	except:
		print(target1, target2)

def read_fasta(path):
	names=[]
	seqs = []
	with open(path, 'r') as f:
		for line in f:
			if line.startswith('>'):
				names.append(line.strip().replace('>', ''))
			else:
				seqs.append(line.strip())
	return zip(names, seqs)

def check_and_create_folder(db_name):
	if not os.path.exists(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name)):
		os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name))

def getamino_uniprot(proteinID):
    r = requests.get(f'https://www.uniprot.org/uniprot/{proteinID}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return proteinID, aminoseq


######################################## START MAIN #########################################
#############################################################################################



def main():
	level= logging.INFO
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)
	parser = argparse.ArgumentParser()
	parser.add_argument("dbPath", help="Path to the database interaction lits",
						default='/home/margaret/data/jfuente/DTI/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv',
						type=str)
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
						help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	args = parser.parse_args()
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

	DB_PATH = args.dbPath
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = get_DB_name(DB_PATH)
	targets = read_and_extract_targets(DB_PATH)
	targets = list(set(targets))
	logging.info('targets: {len(targets)}')
	target_seqs = list(map(getamino_uniprot, tqdm(targets)))
	
	# remove the targets without sequences
	logging.info('Downloading and removing targets without sequences')
	target_seqs = list(filter(lambda x: x[1] != '', target_seqs))

	file_path = os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name, 'Targets_AA_sequences.tsv')
	if not os.path.exists(file_path):
		with open(file_path, 'w') as f:
			for target, seq in target_seqs:
				_ = f.write('>'+target+'\n'+seq+'\n')

	# get the SW scores
	global PATH
	PATH = create_remove_tmp_folder(os.path.join('/tmp/SmithWaterman' , db_name))
	print(PATH)
	all_SmithWaterman = []
	for pair1 in tqdm(target_seqs):
		tmp = []
		if not pair1[1]:
			logging.info(f'No sequence for {pair1[0]}')
			continue
		tmp.extend(repeat(pair1, len(target_seqs)))
		with mp.Pool(processes=mp.cpu_count()-5) as pool:
			results = pool.starmap(get_SW_score, zip(tmp, target_seqs))
		all_SmithWaterman.append(results)

	targets = [ target for target, _ in target_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)
	logging.info('Saving the array')
	check_and_create_folder(db_name)
	file_path = os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name, 'Drugs_SmithWaterman_scores.tsv')
	logging.info('Raw scores saved to: {file_path}')
	SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(PATH)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name, 'Drugs_SmithWaterman_scores_MinMax.tsv')
	logging.info('Normalized scores saved to: {file_path}')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################