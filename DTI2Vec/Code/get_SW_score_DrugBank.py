import argparse
import logging
import requests
import os, sys, uuid
import pandas as pd
import subprocess as sp
import multiprocessing as mp
from itertools import repeat
from re import search
from tqdm import tqdm
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR']
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

def read_and_extract_targets(path):
	dti_file = os.path.join(path, 'DrugBank_DTIs.tsv')
	targets = []
	try:
		with open(dti_file, 'r') as f:
			_ = next(f)
			targets = f.readlines()
			targets = [drug.strip().split('\t')[1] for drug in tqdm(targets, desc='Reading DrugBank DB')]
			return targets
	except FileNotFoundError:
		logging.error('File not found: {}'.format(dti_file))
		sys.exit('Please create the file by running the "get_DTI_from_DrugBank.py" script')

def getamino_KEGG(protein):
    r = requests.get(f'http://rest.kegg.jp/get/{protein}/aaseq')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

def retrieve_sequences(ID):
	"""
	This function retrieves the AA sequence from uniprot.
	"""
	r = requests.get(f'https://www.uniprot.org/uniprot/{ID}.fasta')
	if r.status_code == 200:
		aminoseq = ''.join(r.text.split('\n')[1:])
	else:
		logging.error(f'{ID} not found')
		return (ID, None)
	return (ID, aminoseq)

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
		# logging.debug(f'File {fl_name} already exists')
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
	global already_written_fastas
	target1, seq1 = pair1
	target2, seq2 = pair2
	
	if target1 in already_written_fastas:
		fasta1 = already_written_fastas.get(target1, None)
	else:
		fasta1 = os.path.join(PATH, target1.replace(':', '_')+'.fasta')
		fasta1 = write_fasta(PATH, target1, seq1)
		already_written_fastas[target1] = fasta1
	
	if target2 in already_written_fastas:
		fasta2 = already_written_fastas.get(target2, None)
	else:
		fasta2 = os.path.join(PATH, target2.replace(':', '_')+'.fasta')
		fasta2 = write_fasta(PATH, target2, seq2)
		already_written_fastas[target2] = fasta2
	result_ID = str(uuid.uuid4())
	result_file = os.path.join(PATH, result_ID+'_results.txt')
	args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
			'-asequence', fasta1 , '-bsequence', fasta2, 
			'-gapopen', '10.0', '-gapext', '0.5', 
			'-outfile', result_file]
	try:
		_ = sp.check_call(args, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		score = extract_score(result_file)
		return score
	except:
		logging.warning(f'Not able to compute SW score for : {target1}, {target2}')

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

######################################## START MAIN #########################################
#############################################################################################



def main():
	level= logging.INFO
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)
	parser = argparse.ArgumentParser()
	parser.add_argument("dbPath", help="Path to the database interaction lits",
						default='/home/margaret/data/jfuente/DTI/Data/Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt',
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

	# sanity check for the DB
	# DB_PATH = '/home/margaret/data/jfuente/DTI/Data/DrugBank/'
	DB_PATH = args.dbPath
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = get_DB_name(DB_PATH)

	# read the DB and keep the targets
	targets = read_and_extract_targets(DB_PATH)
	targets = list(set(targets))

	# check if already donwloaded sequences
	folder_path = os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/Data', db_name)
	file_path = os.path.join(folder_path, 'Targets_AA_sequences.tsv')
	if os.path.isfile(file_path):
		logging.info(f'File {file_path} already exists')
		targets_seqs = list(read_fasta(file_path))
	else:
		# get the AA sequences from UniProt
		targets_seqs = []
		for target in tqdm(targets, desc='Getting AA sequences from UniProt'):
			targets_seqs.append(retrieve_sequences(target))

		# remove entries without sequences
		targets_seqs = list(filter(lambda entry: entry[1] != '', targets_seqs))
		targets_seqs = [target for target in targets_seqs if target[1]]
		# write the sequences to fasta
		if not os.path.exists(folder_path):
			os.makedirs(folder_path)
		with open(file_path, 'w') as f:
			for target, seq in targets_seqs:
				_ = f.write('>'+target+'\n'+seq+'\n')

	# get the SW scores
	PATH = create_remove_tmp_folder(os.path.join('/tmp/SmithWaterman' , db_name))
	print(PATH)
	already_written_fastas = {}
	all_SmithWaterman = []
	for pair1 in tqdm(targets_seqs):
		tmp = []
		if not pair1[1]:
			logging.info(f'No sequence for {pair1[0]}')
			continue
		tmp.extend(repeat(pair1, len(targets_seqs)))
		with mp.Pool(processes=mp.cpu_count()-5) as pool:
			results = pool.starmap(get_SW_score, zip(tmp, targets_seqs))
		all_SmithWaterman.append(results)


	targets = [ target for target, _ in targets_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)
	logging.info('Saving the array')
	check_and_create_folder(db_name)
	file_path = os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name, 'Drugs_SmithWaterman_scores.tsv')
	SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(PATH)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/', db_name, 'Drugs_SmithWaterman_scores_MinMax.tsv')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################