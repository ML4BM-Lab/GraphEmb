import os 
import argparse
import logging
import helper_functions_DTI2Vec as hf



######################################## START MAIN #########################################
#############################################################################################
def main():
	level= logging.INFO
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)
	parser = argparse.ArgumentParser()
	parser.add_argument("dbPath", help="Path to the database interaction lits",
						default='./../../DB/Data/DrugBank/DrugBank_DTIs.tsv',
						type=str)
	parser.add_argument("-v", "--verbose", dest="verbosity", default=3,
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
	level= log_levels[int(args.verbosity)]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)
	# we need to create a adjacency matrix with the drugs on the columns and the proteins on the rows
	DB_PATH = args.dbPath
	# DB_PATH = './../../DB/Data/DrugBank/DrugBank_DTIs.tsv'
	###################
	# sanity check for the DB
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
	edges = hf.get_dti(DB_PATH)
	admat = hf.get_admat_from_dti(edges)
	logging.debug(f"Output files at : {os.path.join('./../Data',db_name, db_name +'_admat.tsv')}")
	admat.to_csv(os.path.join('./../Data',db_name ,db_name +'_admat.tsv'), sep='\t')
	hf.write_edges(edges, os.path.join('./../Data',db_name , db_name +'_dti.tsv'))

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################

