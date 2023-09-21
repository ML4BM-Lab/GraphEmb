import os 
import argparse
import logging
from tqdm import tqdm
import multiprocessing as mp
from itertools import repeat
from sklearn.preprocessing import MinMaxScaler
from shutil import rmtree
import pandas as pd
import Protein_kernels_helpers as hf


def YAMANASHI_SW(subdataset='E',model_name='DDR'):
	
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	DB_PATH = './DB/Data/Yamanashi_et_al_GoldStandard/'+subdataset+'/interactions/'+subdataset.lower()+'_admat_dgc_mat_2_line.txt'

	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	# read the DB and keep the targets
	targets = hf.read_and_extract_targets(DB_PATH)
	targets = list(set(targets))

	file_path = os.path.join('.',model_name,'Data/Yamanashi_et_al_GoldStandard',subdataset,subdataset+'_Targets_AA_sequences.tsv')

	if os.path.isfile(file_path):
		logging.info(f'Reading AA sequences from: {file_path}')
		targets_seqs = hf.read_AA_sequences(file_path)
	else:
		# get the AA sequences from KEGG
		targets_seqs = []
		for target in tqdm(targets, desc='Getting AA sequences from KEGG'):
			targets_seqs.append(hf.replace_and_get_AA(target))
		# remove entries without sequences
		targets_seqs = list(filter(lambda entry: entry[1] != '', targets_seqs))
		# write the sequences to fasta
		if not os.path.exists(file_path):
			with open(file_path, 'w') as f:
				for target, seq in targets_seqs:
					_ = f.write('>'+target+'\n'+seq+'\n')

	# get the SW scores
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
	
	targets = [ target for target, _ in targets_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)
	logging.info('Saving the array')
	hf.check_and_create_folder(db_name,model_name)
	file_path = os.path.join('.',model_name,'Data', db_name, subdataset+'_prot_SmithWaterman_scores.tsv')
	#SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(tmp_path)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('.',model_name,'Data', db_name, subdataset+'_prot_SmithWaterman_scores_MinMax.tsv')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')

def BIOSNAP_SW(model_name='DDR', subdataset='BIOSNAP'):
    
	fmt = '[%(levelname)s] %(message)s]'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	DB_PATH = './DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
	#DB_PATH = '/mnt/md0/data/jfuente/DTI/Input4Models/DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'

	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	targets = hf.read_and_extract_targets(DB_PATH)
	targets = list(set(targets))
	logging.info(f'targets: {len(targets)}')
	targets_seqs = list(map(hf.retrieve_sequences_from_UniProt, tqdm(targets)))
	
	# remove the targets without sequences
	logging.info('Downloading and removing targets without sequences')
	targets_seqs = list(filter(lambda x: x[1] != '' and x[1] != None, targets_seqs))

	file_path = os.path.join('.',model_name,'Data/BIOSNAP',subdataset+'_Targets_AA_sequences.tsv')
	#file_path = os.path.join('/mnt/md0/data/jfuente/DTI/Input4Models',model_name,'Data/BIOSNAP',subdataset+'_Targets_AA_sequences.tsv')

	if not os.path.exists(file_path):
		with open(file_path, 'w') as f:
			for target, seq in tqdm(targets_seqs,desc='Writting files'):
				_ = f.write('>'+target+'\n'+seq+'\n')

	# get the SW scores
	tmp_path  = hf.create_remove_tmp_folder(os.path.join('/media/scratch_ssd/tmp/' , db_name))
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


	targets = [ target for target, _ in targets_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)
	logging.info('Saving the array')
	hf.check_and_create_folder(db_name,model_name)
	file_path = os.path.join('.', model_name,'Data', subdataset, subdataset+'_prot_SmithWaterman_scores.tsv')
	logging.info('Raw scores saved to: {file_path}')
	#SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(tmp_path)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('.', model_name, 'Data', db_name, subdataset+'_prot_SmithWaterman_scores_MinMax.tsv')
	logging.info('Normalized scores saved to: {file_path}')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')


def DRUGBANK_SW(model_name = 'DDR', subdataset = 'DrugBank'):

	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	# sanity check for the DB
	DB_PATH = './DB/Data/DrugBank/DrugBank_DTIs.tsv'

	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	# read the DB and keep the targets
	targets = hf.read_and_extract_targets(DB_PATH)
	targets = list(set(targets))

	# check if already donwloaded sequences
	file_path = os.path.join('.',model_name,'Data/DrugBank',subdataset+'_Targets_AA_sequences.tsv')

	if os.path.isfile(file_path):
		logging.info(f'File {file_path} already exists')
		targets_seqs = list(hf.read_fasta(file_path))
	else:
		# get the AA sequences from UniProt
		targets_seqs = []
		for target in tqdm(targets, desc='Getting AA sequences from UniProt'):
			targets_seqs.append(hf.retrieve_sequences_from_UniProt(target))

		# remove entries without sequences
		targets_seqs = list(filter(lambda entry: entry[1] != '', targets_seqs))
		targets_seqs = [target for target in targets_seqs if target[1]]
		# write the sequences to fasta
		if not os.path.exists(file_path):
			with open(file_path, 'w') as f:
				for target, seq in targets_seqs:
					_ = f.write('>'+target+'\n'+seq+'\n')

	# get the SW scores
	tmp_path  = hf.create_remove_tmp_folder(os.path.join('/media/scratch_ssd/tmp/' , db_name))
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

	targets = [ target for target, _ in targets_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)
	logging.info('Saving the array')
	hf.check_and_create_folder(db_name,model_name)
	file_path = os.path.join('.', model_name,'Data', subdataset, subdataset+'_prot_SmithWaterman_scores.tsv')
	#SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(tmp_path)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('.', model_name, 'Data', db_name, subdataset+'_prot_SmithWaterman_scores_MinMax.tsv')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')

def DAVIS_SW(model_name = 'DDR', subdataset = 'Davis_et_al'):

	fmt = '[%(levelname)s] %(message)s'
	
	logging.basicConfig(format=fmt, level=logging.DEBUG)
	DB_PATH =  './DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al.tsv'

	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	
	targets_seqs = list(set(hf.get_seqs_DAVIS(DB_PATH)))
	logging.info(f'{len(targets_seqs)} targets found')

	file_path = os.path.join('.',model_name,'Data/Davis_et_al/Davis_et_al_Targets_AA_sequences.tsv')

	if not os.path.isfile(file_path):
		with open(file_path, 'w') as f:
				for target, seq in targets_seqs:
					_ = f.write('>'+target+'\n'+seq+'\n')
	
	# get the SW scores
	tmp_path  = hf.create_remove_tmp_folder(os.path.join('/media/scratch_ssd/tmp/' , db_name))
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

	targets = [ target for target, _ in targets_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)
	logging.info('Saving the array')
	hf.check_and_create_folder(db_name,model_name)
	file_path = os.path.join('.', model_name,'Data', db_name, subdataset+'_prot_SmithWaterman_scores.tsv')
	#SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(tmp_path)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('.', model_name, 'Data', db_name, subdataset+'_prot_SmithWaterman_scores_MinMax.tsv')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')

def BINDINGDB_SW(model_name = 'DDR', dataset = 'BindingDB'):

	fmt = '[%(levelname)s] %(message)s]'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	DB_PATH =  './DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
	#DB_PATH = args.dbPath
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	targets_seqs = list(set(hf.get_seqs_BindingDB(DB_PATH)))
	file_path = os.path.join('.',model_name,'Data/BindingDB/BindingDB_Targets_AA_sequences.tsv')

	if not os.path.exists(file_path):
		with open(file_path, 'w') as f:
			for target, seq in targets_seqs:
				_ = f.write('>'+target+'\n'+seq+'\n')
	# get the SW scores
	tmp_path  = hf.create_remove_tmp_folder(os.path.join('/media/scratch_ssd/tmp/' , db_name))
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


	targets = [ target for target, _ in targets_seqs ]
	SmithWaterman_arr = pd.DataFrame(all_SmithWaterman, columns=targets, index=targets)
	logging.info('Saving the array')
	hf.check_and_create_folder(db_name,model_name)
	file_path = os.path.join('.', model_name,'Data', db_name, dataset+'_prot_SmithWaterman_scores.tsv')
	logging.info('Raw scores saved to: {file_path}')
	#SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(tmp_path)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	#file_path = os.path.join('./../Data', db_name, '_prot_SmithWaterman_scores_MinMax.tsv')
	file_path = os.path.join('.', model_name, 'Data', db_name, dataset+'_prot_SmithWaterman_scores_MinMax.tsv')
	logging.info('Normalized scores saved to: {file_path}')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')

