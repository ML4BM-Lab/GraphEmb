import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from tqdm import tqdm
import logging
import os
import helper_functions_dtinet as hf
import argparse
import helper_splits_eegdti as splits
from tqdm.contrib.itertools import product

def main():
	'''
	generate one to one index foldere & files
	considering splits and subsampling 
	'''
	parser = argparse.ArgumentParser() 
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=4,
					help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	parser.add_argument("-dbPath","--dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)
	parser.add_argument("-split_type", "--split_type", help="Select the type of split ['Sp', 'Sd', 'St'] to generate oneTooneIndex folder", type=str)
	parser.add_argument("-subsampling", help="Flag for subsampling True", action='store_true')

	args = parser.parse_args()

	# -) log info; 
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

	#######
	### log output detals
	logging.info(
		'''
		This script generates splits
		'''
		)
	# OUTPUT DIRECTORY
	# sanity check
	DB_PATH = args.dbPath
	SPLIT_TYPE = args.split_type
	if SPLIT_TYPE not in  ['Sp', 'Sd', 'St']:
		raise NameError('Need to specify a valid split type')

	SUBSAMPLING_TYPE = args.subsampling
	#if DB_PATH not in  [True, False] :
	#	raise NameError('Need to specify a subsampling typre as True or False')

	logging.info(f'Working in output folder for: {DB_PATH}')
	logging.info(f'Creating Split: {SPLIT_TYPE}')
	logging.info(f'Subsampling: {SUBSAMPLING_TYPE}')

	db_name = hf.get_DB_name(DB_PATH)
	hf.check_and_create_folder(db_name)
	# Create relative output path
	wdir = os.path.join('../Data', db_name)
	# wdir = '../Data/DrugBank'

	##########################################
	# create oneTooneIndex folder if it does not exists
	# type Sp first
	#if SPLIT_TYPE == 'Sp':
	if SUBSAMPLING_TYPE:
		sub = 'subsampling'
	else:
		sub = 'nosubsampling'
	path_folder = os.path.join(wdir, f'oneTooneIndex_{SPLIT_TYPE}_{sub}')
	if not os.path.exists(path_folder):
		os.makedirs(path_folder)

	fpath = os.path.join(os.getcwd(), wdir, f'final_dtis_{DB_PATH}.tsv')
	DTIs = pd.read_csv(fpath,sep='\t',header=None)
	DTIs.columns = ['Drug','Protein']
	logging.debug(DTIs.head())
	sp_splits = splits.generate_splits(DTIs, mode= SPLIT_TYPE, subsampling=SUBSAMPLING_TYPE, foldnum=10)
	logging.debug(sp_splits[0][0][0][0])
	# this also changes (drug, protein) to (protein, drug)
	sp_splits = splits.get_id2idx_protdrug(wdir, sp_splits)
	logging.debug(sp_splits[0][0][0][0])
	nseed, nfold = 0, 0 
	for nseed, nfold in product(range(len(sp_splits)), range(len(sp_splits[nseed]))):
		np.savetxt(os.path.join(path_folder, f'train_index_pos(0,1){nseed}_{nfold}.txt'), sp_splits[nseed][nfold][0], fmt='%i', delimiter=" ")
		np.savetxt(os.path.join(path_folder, f'train_index_neg(0,1){nseed}_{nfold}.txt'), sp_splits[nseed][nfold][1], fmt='%i', delimiter=" ")
		np.savetxt(os.path.join(path_folder, f'test_index_pos(0,1){nseed}_{nfold}.txt'), sp_splits[nseed][nfold][2], fmt='%i', delimiter=" ")
		np.savetxt(os.path.join(path_folder, f'test_index_neg(0,1){nseed}_{nfold}.txt'), sp_splits[nseed][nfold][3], fmt='%i', delimiter=" ")
	

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################


