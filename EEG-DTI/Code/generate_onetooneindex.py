import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from tqdm import tqdm
import logging
import os
import helper_functions_dtinet as hf
import argparse



def main():
	'''
	DrugBank
	'''
	parser = argparse.ArgumentParser() 
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=4,
					help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	parser.add_argument('-json', help="If selected, outputs a dictionary in json", action="store_true")
	parser.add_argument("dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)
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
		This script generates ....
		'''
		)
	# OUTPUT DIRECTORY
	# sanity check
	DB_PATH = args.dbPath
	logging.info(f'Working in output folder for: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	hf.check_and_create_folder(db_name)
	# Create relative output path
	wdir = os.path.join('../Data', db_name)
	# wdir = '../Data/DrugBank'
	##########################################
	# create oneTooneIndex folder if it does not exists
	folder = 'oneTooneIndex'
	path_folder = os.path.join(wdir, folder)
	if not os.path.exists(path_folder):
		os.makedirs(path_folder)
	# read the adjacency matrix from sevenets
	# protein drug
	admat = pd.read_csv(os.path.join(wdir, 'sevenNets/mat_protein_drug.txt'), sep=' ', header=None)
	logging.debug(f'admat shape: {admat.shape}')
	# select edges
	rows = admat.index.tolist()
	columns = admat.columns.tolist()
	edges_all_ = []
	edges_all_false_ = []
	for row in tqdm(rows):
		for column in columns:
			element = admat.loc[row][column]
			if element == 1:
				onetoone = (rows.index(row), columns.index(column))
				#print(onetoone)
				edges_all_.append(onetoone)
			elif element == 0:
				onetoone = (rows.index(row), columns.index(column))
				edges_all_false_.append(onetoone)
			else:
				logging.info('weird edges in binary matrix')
	edges_all = pd.DataFrame(edges_all_)
	logging.info(f'# edges postives: {edges_all.shape[0]} ')
	edges_all_false = pd.DataFrame(edges_all_false_)
	logging.info(f'# edges false: {edges_all_false.shape[0]} ')
	# edges false need to be cut with the same shape as edges true, randomly
	# check the section where we load edges_false_protein_drug_01 to see this
	edges_all_false = edges_all_false.sample(n=edges_all.shape[0], random_state=123)
	logging.info(f'# edges false now with sample (random state 123) : {edges_all_false.shape[0]} ')
	# this is how they "balance" the dataset, ramdomly chosing negatives
	# save
	logging.debug('writting all edges files')
	edges_all.to_csv(os.path.join(path_folder, 'edges_all_(0, 1).txt'), sep= " ", header=None, index=None)
	edges_all_false.to_csv(os.path.join(path_folder, 'edges_all_false(0, 1).txt'), sep= " ", header=None, index=None)
	# 
	all_edge_idx = list(range(edges_all.shape[0]))
	np.random.seed(123)
	np.random.shuffle(all_edge_idx) 
	kf = KFold(n_splits=10, shuffle=True, random_state=123)
	# Then we use the 10-flod cross-validation method to evaluate our mode
	# we first generate the train set and test set.
	i =0
	for pos_train_idx, pos_test_idx in kf.split(all_edge_idx):
		# pos_train_idx is train set, pos_test_idx is test set
		###### TRAIN  INDEX ( 0, 1 )
		train_edges = edges_all.iloc[pos_train_idx]  # These are  'train_index_(0, 1)0.txt'
		train_edges.to_csv(os.path.join(path_folder, f'train_index_(0, 1){i}.txt'), sep= " ", header=None, index=None)
		#
		###### TEST  INDEX ( 0, 1 )
		test_edges = edges_all.iloc[pos_test_idx] 
		logging.debug(f'test edge shape: {test_edges.shape}')
		test_edges.to_csv(os.path.join(path_folder, f'test_index_(0, 1){i}.txt'), sep= " ", header=None, index=None)
		i+=1
	# 
	### negative test splits?
	# we dont know exactly how they did
	# so we repeat the process with the negative matrix a
	# and only sabe those for test (only txt that the model takes)
	all_edge_idx_false = list(range(edges_all_false.shape[0]))
	np.random.seed(123)
	np.random.shuffle(all_edge_idx_false) # ? 
	#all_edge_idx_false
	# define kf again but actually not needed
	# gonna be the same as the random state is 123
	kf = KFold(n_splits=10, shuffle=True, random_state=123)
	i =0
	for pos_neg_train_idx, pos_neg_test_idx in kf.split(all_edge_idx_false):
		# They dont use the train neg edges 
		#train_neg_edges = edges_all_false.iloc[pos_neg_train_idx]  # The
		#print(train_neg_edges.shape)
		test_neg_edges = edges_all_false.iloc[pos_neg_test_idx]  # The
		logging.debug(f'test neg edge shape: {test_neg_edges.shape}')
		test_neg_edges.to_csv(os.path.join(path_folder, f'index_test_false(0, 1){i}.txt'), sep= " ", header=None, index=None)
		i+=1


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################


