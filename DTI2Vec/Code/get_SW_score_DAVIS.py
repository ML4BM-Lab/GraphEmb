import argparse
from asyncio.log import logger
import logging
import os
import pandas as pd
from tqdm import tqdm
from itertools import repeat
import multiprocessing as mp
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler
import helper_functions_DTI2Vec as hf

######################################## START MAIN #########################################
#############################################################################################


def main():
	level= logging.INFO
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)
	parser = argparse.ArgumentParser()
	parser.add_argument("dbPath", help="Path to the database interaction lits",
						default='/home/margaret/data/jfuente/DTI/Data/Davis_et_al/Davis_30.tsv',
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
	fmt = '[%(levelname)s] %(message)s]'
	logging.basicConfig(format=fmt, level=level)
	# DB_PATH =  './../../DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al.tsv'
	DB_PATH = args.dbPath
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	
	targets_seqs = list(set(hf.get_seqs_DAVIS(DB_PATH)))
	logging.info(f'{len(targets_seqs)} targets found')
	file_path = os.path.join('./../Data/', db_name, 'Targets_AA_sequences.tsv')
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
	hf.check_and_create_folder(db_name)
	file_path = os.path.join('./../Data', db_name, 'Proteins_SmithWaterman_scores.tsv')
	SmithWaterman_arr.to_csv(file_path, sep='\t')
	rmtree(tmp_path)
	zscore_SmithWaterman_arr = pd.DataFrame(MinMaxScaler().fit_transform(SmithWaterman_arr),columns=targets,index=targets)
	file_path = os.path.join('./../Data', db_name, 'Proteins_SmithWaterman_scores_MinMax.tsv')
	zscore_SmithWaterman_arr.to_csv(file_path, sep='\t')




#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################