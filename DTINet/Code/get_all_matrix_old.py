import os, sys
import pandas as pd
import logging
import argparse
from re import search

'''
This script should return al mat_A_B.txt
given all coordinate files (tsv) as input
and XXXXX
'''


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank', 'E', 'GPCR', 'IC', 'NR']
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


def check_and_create_folder(db_name):
	if not os.path.exists(os.path.join('../Data', db_name)):
		os.mkdir(os.path.join('../Data', db_name))


wdir = '../Data/BindingDB'
list_of_files = ['coordinates_drug_disease.tsv', 'coordinates_drug_se.tsv', 
                'coordinates_protein_disease.tsv', 'coordinates_drug_drug.tsv', 
                'coordinates_PPI.tsv']
# all coordinate files should have the first column the element that appear as row
# all coordinate files should have as the 2nd colum the element that appear as columns
drug_dis = pd.read_csv(os.path.join(wdir, list_of_files[0]), sep='\t')
drug_se = pd.read_csv(os.path.join(wdir, list_of_files[1]), sep='\t')
PPI = pd.read_csv(os.path.join(wdir, list_of_files[2]), sep='\t', names=['P1', 'P2'])
drug_drug = pd.read_csv(os.path.join(wdir, list_of_files[3]), sep='\t', names=['D1', 'D2'])
DTI = pd.read_csv(os.path.join(wdir, list_of_files[4]), sep='\t', names=['Drug', 'Protein'])
prot_dis = pd.read_csv(os.path.join(wdir, list_of_files[5]), sep='\t', usecols=['UniprotID', 'DiseaseID'])




# mal: coordinates_drug_se => change & make DrugBank_ID//se ---> cambiando


# selection for drug nodes (for index/row)


# selection for protein (for index/row)

######################################## START MAIN #########################################
#############################################################################################
def main():
    '''
    Getting all matrix  
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
    parser.add_argument('-json', help="If selected, outputs a dictionary in json", action="store_true")
    parser.add_argument("dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)
   #
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
    #
    #######
    ###  output details
    logging.info(
        '''
        This script needs:
            - From SIDER: meddra_all_se.tsv, 
            - DrugBank xml file


        This script generates the following files:
            - se.txt
            - coordinates_drug_se.tsv
        
        DTINet uses from CTD the processed files:
            - mat_drug_se.txt
        '''
        )
    # WORKING DIRECTORY
    DB_PATH = args.dbPath
    db_name = get_DB_name(DB_PATH)
    check_and_create_folder(db_name)
    # Create relative output path ? outputpath --> working directory
    wdir = os.path.join('../Data', db_name)
    ##########################################
    # READING ALL TSV FILES
    
