# import sys
# sys.path.append("/home/margaret/data/jfuente/DTI/Code/CommonScripts/")
import os, uuid
import numpy as np
from GetDrugsAndGenes import get_SIMMCOMP_score, getamino_KEGG
import multiprocessing as mp
from itertools import repeat
from tqdm import tqdm
import subprocess as sp
import logging

level= logging.INFO
fmt = '[%(levelname)s] %(asctimes)s - %(message)s'
logging.basicConfig(format=fmt, level=level)

def get_all_DT(DB_PATH):
	"""
	This function returns all Drugs and Targets from a DB wit 'Drug\tTarget'.
	"""
	all_drugs = []
	all_targets = []
	# list all DB 
	all_DB = os.listdir(DB_PATH)
	for db in all_DB:
		db_file = os.path.join(DB_PATH, db, 'interactions')
		file_name = os.listdir(db_file)[0]
		with open(os.path.join(db_file, file_name), 'r') as f:
			for line in f:
				drug, target = line.split('\t')
				all_drugs.append(drug)
				all_targets.append(target.strip())
	return all_drugs, all_targets


DB_PATH = '/home/margaret/data/jfuente/DTI/Data/Yamanashi_et_al_GoldStandard'


drugs, targets = get_all_DT(DB_PATH)
# Keep the unique drugs and targets
drugs   = list(set(drugs))
targets = list(set(targets))


'''
===========================================================================================
								SIMCOMP calculation
===========================================================================================
'''
# Compute the score for all pairs of drugs and create the array
# ======   using multiporcess
# Pool = mp.Pool(5)
# all_SIMCOMP_scores=[]
# for query1 in tqdm(drugs[0:10], desc='Retrieving scores from KEGG'):
# 	for query2 in drugs[0:10]:
# 		tmp =[]
# 		tmp.extend(repeat(query1, len(drugs[0:10])))
# 		results = Pool.starmap(get_SIMMCOMP_score, zip(tmp, drugs[0:10]))
# 		results = list(map(lambda x: float(x.split('\t')[-1].strip()), results))
# 	all_SIMCOMP_scores.append(results)

# SIMCOMP_arr = np.asarray(all_SIMCOMP_scores)

# ======   sequential
all_SIMCOMP_scores=[]
for query1 in tqdm(drugs, desc='Retrieving scores from KEGG'):
	for query2 in tqdm(drugs):
		results = get_SIMMCOMP_score(query1, query2)
		results = results.split('\t')[-1].strip()
	all_SIMCOMP_scores.append(results)

SIMCOMP_arr = np.asarray(all_SIMCOMP_scores)

# Get the AA sequences of the targets
# WARNING: the protein identifier needs to be in the format 'hsa:'+number
# targets_seqs = list(map(lambda identifier: (identifier, getamino_KEGG(identifier.replace('hsa', 'hsa:'))), targets))




'''
===========================================================================================
								SmithWaterman Scores
===========================================================================================
'''

def replace_and_get_AA(ID):
	"""
	This function replaces the 'hsa' in the identifier with 'hsa:'.
	"""
	ID = ID.replace('hsa', 'hsa:')
	return (ID, getamino_KEGG(ID))

# targets_seqs = Pool.map(replace_and_get_AA, targets)

targets_seqs = []
for target in tqdm(targets):
	targets_seqs.append(replace_and_get_AA(target))

with open('/home/margaret/data/jfuente/DTI/InputData/DTI2Vec/Yamanashi_et_al_GoldStandard/Target_sequences.fasta', 'w') as f:
	for target, seq in targets_seqs:
		_ = f.write('>'+target+'\n'+seq+'\n')

###-------------------------- USE BLAST --------------------------------###
USE_BLAST = False
if USE_BLAST:
	# save the files to fasta
	# path_2_fastas = '/home/margaret/data/jfuente/DTI/Data/Yamanashi_et_al_GoldStandard/targets_seqs.fasta'
	path_2_fastas = '/tmp/targets_seqs.fasta'
	# with open(path_2_fastas, 'w') as f:
	with open(path_2_fastas, 'w') as f:
		for target, seq in targets_seqs:
			_ = f.write('>'+target+'\n'+seq+'\n')
	data_set_name = 'Yamanashi_et_al_GoldStandard'
	# Create the custom DB
	path_2_makeDB = '/home/margaret/data/gserranos/BLAST_StandAlone/ncbi-blast-2.12.0+/bin/makeblastdb'
	args = [path_2_makeDB, '-dbtype', 'prot', '-in', path_2_fastas, '-input_type', 'fasta', '-title', data_set_name, '-out', data_set_name]
	sp.check_call(args)
	path_2_blastp = '/home/margaret/data/gserranos/BLAST_StandAlone/ncbi-blast-2.12.0+/bin/blastp'
	args = [path_2_blastp, '-db', data_set_name, '-query', path_2_fastas, '-outfmt', '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore score', '-out', data_set_name+'.blast_results']
	sp.check_call(args)


###-------------------------- USE SMITH WATERMAN --------------------------------###

def extract_score(results_file):
	with open(results_file, 'r') as f:
		for line in f:
			if not line.startswith('# Score:'):
				continue
			else:
				return float(line.split()[-1])

def create_remove_tmp_folder(path):
	if os.path.exists(path):
		logging.info('Removing tmp folder: {}'.format(path))
		os.system('rm -r '+path)
	else:
		logging.info('Creating tmp folder: {}'.format(path))
		os.mkdir(path)
		return path

def write_fasta(path, target, seq):
	fl_name = os.path.join(path, target.replace(':', '_')+'.fasta')
	with open(fl_name, 'w') as f:
		_ = f.write('>'+target+'\n'+seq+'\n')
	return fl_name

path = create_remove_tmp_folder('/tmp/SmithWaterman')


path_2_water = '/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water'
all_SmithWaterman = []
for target1, seq1 in tqdm(targets_seqs):
	if not seq1:
		continue
	path1 = write_fasta(path, target1, seq1)
	tmp_results=[]
	for target2, seq2 in tqdm(targets_seqs):
		if not seq2:
			continue
		path2 = write_fasta(path, target2, seq2)
		result_file = os.path.join(path, 'results.txt')
		args = [path_2_water, '-asequence', path1 , 
				'-bsequence', path2, '-gapopen', '10.0', 
				'-gapext', '0.5', '-outfile', result_file]
		_ = sp.check_call(args, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		tmp_results.append(extract_score(result_file))
	all_SmithWaterman.append(tmp_results)




all_SmithWaterman = np.asarray(all_SmithWaterman)


def get_SW_score(pair1, pair2):
	target1, seq1 = pair1
	target2, seq2 = pair2
	if not seq2:
		return 0
	fasta1 = os.path.join(path, target1.replace(':', '_')+'.fasta')
	if not os.path.exists(fasta1):
		fasta1 = write_fasta(path, target1, seq1)
	fasta2 = write_fasta(path, target2, seq2)
	result_ID = str(uuid.uuid4())
	result_file = os.path.join(path, result_ID+'_results.txt')
	args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
			'-asequence', fasta1 , '-bsequence', fasta2, 
			'-gapopen', '10.0', '-gapext', '0.5', 
			'-outfile', result_file]
	try:
		_ = sp.check_call(args, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		return extract_score(result_file)
	except:
		print(target1, target2)




path_2_water = '/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water'
all_SmithWaterman = []
for pair1 in tqdm(targets_seqs):
	tmp = []
	if not pair1[1]:
		continue
	tmp.extend(repeat(pair1, len(targets_seqs)))
	with mp.Pool(processes=mp.cpu_count()) as pool:
		results = pool.starmap(get_SW_score, zip(tmp, targets_seqs))
	all_SmithWaterman.append(results)
