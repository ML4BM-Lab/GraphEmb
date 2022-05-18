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

def BIOSNAP_SW():
    pass

def DRUGBANK_SW():
    pass

def DAVIS_SW():
    pass

def BINDINGDB_SW():
    pass

