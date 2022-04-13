import requests
import os, sys, uuid, re
import pandas as pd
import subprocess as sp
import multiprocessing as mp
from itertools import repeat
from re import search
from tqdm import tqdm
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler
import logging

def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'DrugBank' ,'E', 'GPCR', 'IC', 'NR']
	for db in DB_NAMES:
		if search(db, path):
			logging.info(f'Database: {db}')
			if db in ['E', 'GPCR', 'IC', 'NR']:
				db = os.path.join('Yamanashi_et_al_GoldStandard', db)
				if not os.path.exists(os.path.join('./../Data', db)):
					logging.infor(f'Creating direcotry: ./../Data/{db}')
					os.mkdir(os.path.join('./../Data', db))
				return db
			else:
				if not os.path.exists(os.path.join('./../Data', db)):
					logging.infor(f'Creating direcotry: ./../Data/{db}')
					os.mkdir(os.path.join('./../Data', db))
				return db
	logging.error(f'Database: {db} not found')
	sys.exit('Please provide a valid database')

def read_dtis_Yamanashi(path):
	""""
	Read the dti file and return the dti along with 
	the numerical dictionary sorted by name
	"""
	with open(path, 'r') as f:
		dtis = f.readlines()
	dtis = [entry.strip().split('\t') for entry in dtis]
	all_elements = list(set(element for entry in dtis for element in entry))
	logging.info(f'Number of unique nodes: {len(all_elements)}')
	logging.info(f'Number of unique edges: {len(dtis)}')
	dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
	return dtis, dictionary 

def read_dtis_DrugBank(path):
	""""
	Read the dti file and return the dti along with 
	the numerical dictionary sorted by name
	"""
	with open(path, 'r') as f:
		_ = next(f)
		dtis = f.readlines()
	dtis = [entry.strip().split('\t') for entry in dtis]
	all_elements = list(set(element for entry in dtis for element in entry))
	logging.info(f'Number of unique nodes: {len(all_elements)}')
	logging.info(f'Number of unique edges: {len(dtis)}')
	dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
	return dtis, dictionary 

def read_dtis_BIOSNAP(path):
	""""
	Read the dti file and return the dti along with 
	the numerical dictionary sorted by name
	"""
	with open(path, 'r') as f:
		_ = next(f)
		dtis = f.readlines()
	dtis = [entry.strip().split('\t') for entry in dtis]
	all_elements = list(set(element for entry in dtis for element in entry))
	logging.info(f'Number of unique nodes: {len(all_elements)}')
	logging.info(f'Number of unique edges: {len(dtis)}')
	dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
	return dtis, dictionary 

def read_dtis_DAVIS(path):
	""""
	Read the dti file and return the dti along with 
	the numerical dictionary sorted by name
	"""
	with open(path, 'r') as f:
		_ = next(f)
		dtis = f.readlines()
	dtis = [(entry.split('\t')[1], entry.split('\t')[3])  for entry in dtis]
	all_elements = list(set(element for entry in dtis for element in entry))
	logging.info(f'Number of unique nodes: {len(all_elements)}')
	logging.info(f'Number of unique edges: {len(dtis)}')
	dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
	return dtis, dictionary

def read_dtis_Binding(path):
	""""
	Read the dti file and return the dti along with 
	the numerical dictionary sorted by name
	"""
	with open(path, 'r') as f:
		_ = next(f)
		dtis = f.readlines()
	dtis = [(entry.split('\t')[1], entry.split('\t')[3])  for entry in dtis]
	all_elements = list(set(element for entry in dtis for element in entry))
	logging.info(f'Number of unique nodes: {len(all_elements)}')
	logging.info(f'Number of unique edges: {len(dtis)}')
	dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
	return dtis, dictionary


def write_dtis(dtis, path):
	result_ID = str(uuid.uuid4())
	result_file = os.path.join(path, result_ID+'_coded_dti.txt')
	TAB = '\t'
	NL = '\n'
	dtis = [TAB.join(entry) for entry in dtis]
	with open(result_file, 'w') as outfl:
		outfl.writelines(NL.join(dtis))
	logging.debug(f'Coded DTIs written at {result_file}') 
	return result_file

def write_dict(dictionary, file_path):
	with open(file_path, 'w') as f:
		for key, value in dictionary.items():
			f.write(f'{key}\t{value}\n')
	logging.debug(f'Coded DTIs written at {file_path}') 

def write_embedddings_with_name(embeddings, node_dict):
	with open(embeddings , 'r') as infl :
		stats = next(infl).strip().split()
		embs = [line.strip().split() for line in infl]

	assert len(embs)    == int(stats[0]), 'Number of embeddings does not match'
	assert len(embs[0]) == (int(stats[1]) +1 ), 'Number of dimensions does not match'
	inv_map = {v: k for k, v in node_dict.items()}
	SEP = ' '
	named_nodes = embeddings.replace('.txt', '_with_name.txt')
	with open(named_nodes, 'w') as outfl:
		for emb in embs:
			node = emb[0]
			emb = SEP.join(emb[1:])
			outfl.write(f'{inv_map[node]}{SEP}{emb}\n')
	logging.debug(f'Embeddings written at {named_nodes}')



def create_remove_tmp_folder(path):
	if not os.path.exists(path):
		logging.info('Creating tmp folder: {}'.format(path))
		os.makedirs(path)
		return path
	else: 
		return path

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

def getamino_KEGG(protein):
    r = requests.get(f'http://rest.kegg.jp/get/{protein}/aaseq')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

def replace_and_get_AA(ID):
	"""
	This function replaces the 'hsa' in the identifier with 'hsa:'.
	"""
	ID = ID.replace('hsa', 'hsa:')
	return (ID, getamino_KEGG(ID))

def create_remove_tmp_folder(path):
	if not os.path.exists(path):
		logging.info('Creating tmp folder: {}'.format(path))
		os.makedirs(path)
		return path
	else: 
		return path

def get_yamanashi_subDB(path):
	return re.search(r'(?<=\/)[\w]+', path).group()

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

def get_SW_score(pair1, pair2, tmp_path):
	target1, _ = pair1
	target2, _ = pair2
	fasta1 = os.path.join(tmp_path, target1.replace(':', '_')+'.fasta')
	fasta2 = os.path.join(tmp_path, target2.replace(':', '_')+'.fasta')
	args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
			'-asequence', fasta1 , '-bsequence', fasta2, 
			'-gapopen', '10.0', '-gapext', '0.5', 
			'-stdout']
	try:
		score = sp.Popen(args, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.DEVNULL)
		score = score.communicate(b"\n")
		score = score[0].decode().split('\n')
		score = extract_score(score)
		return score
	except:
		logging.warning(f'Not able to compute SW score for : {target1}, {target2}')

def extract_score(score):
	score = [line for line in score if line.startswith('# Score:')]
	if score:
		return float(score[0].split()[-1])
	else:
		return None

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
	if not os.path.exists(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/', db_name)):
		os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/', db_name))

def read_AA_sequences(path):
	# read a fasta file and return a dictionary with the target as key and the sequence as value
	AA_sequences = {}
	with open(path, 'r') as f:
		for line in f:
			if line.startswith('>'):
				target = line.strip().replace('>', '')
			else:
				AA_sequences[target] = line.strip()
	return list(AA_sequences.items())

def write_all_fastas(fastas, path):
	for header, seq in fastas:
		write_fasta(path, header, seq)
	logging.info('All fastas written')
