import os, sys
import argparse
import logging
import numpy as np
import pandas as pd
import json
import requests
from tqdm import tqdm
from re import search
import xml.etree.ElementTree as ET


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




def get_drug_drug_coordinates(drug_entry, drug_drug_coodinates):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	all_interactions =  drug_entry.findall('.//{http://www.drugbank.ca}drug-interactions')[0]
	for inter in all_interactions: # aqui iteramos en targets
		coordinate = (drugbank_ID, list(inter)[0].text)	
		drug_drug_coodinates.append(coordinate)
	return drug_drug_coodinates


######################################## START MAIN #########################################
#############################################################################################


def main():
    '''
    DrugBank for side information (not DTI) 
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
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
    ###  output details
    logging.info(
        '''
        This script needs:
            - DrugBank xml file

        This script generates the following files:
            - all_drugs.txt
            - all_drug_dic_map.txt (/ dic_drugnames_DBID.json)
            - coordinates_drug_se.tsv
        
        DTINet uses from DrugBank the processed files:
            - mat_drug_drug.txt
        '''
        )

    # OUTPUT DIRECTORY
    DB_PATH = args.dbPath
    db_name = get_DB_name(DB_PATH)
    check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)
    ## SIDER DATA FOLDER
    #data_path = '../../../Data/cross_side_information_DB/SIDER'
    # info 
    #logging.info(f'Processing {data_path[40:]} in {output_path}')

    ### do stuff
    logging.info('Reading DrugBank xml file...')
    tree = ET.parse('../../../Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    root = tree.getroot()

    ########### DRUG NODES ###########

    drug_IDs = []
    drug_names = []
    for drug_entry in tqdm(root):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        name = drug_entry.find('{http://www.drugbank.ca}name').text # incluir aqui name sin check
        drug_names.append(name)
        drug_IDs.append(drugbank_ID)

    # Drug nodes (all_drugs.txt). All nodes in DB (w\o filter) 
    logging.info(f'Writing all_drugs.txt...')
    with open(os.path.join(output_path, 'all_drugs.txt'), 'w') as f:
        #_ = f.write('#DrugBankID\n')
        for item in drug_IDs:
            _ = f.write("%s\n" % item)

    # Write dictionary (drugbank_ID and name)
    logging.info(f'Writing all_drug_dic_map.txt...')
    with open(os.path.join(output_path,'all_drug_dic_map.txt'), 'w') as f:
        for i in range(len(drug_IDs)):
            _ = f.write("%s:%s\n" % (drug_IDs[i],drug_names[i]))

    # json
    if args.json == True:
        dic_drugnames_DBID = dict(zip(drug_IDs, drug_names)) 
        file_name_json = 'dic_drugnames_DBID.json'
        logging.info(f'Writing {file_name_json}...')
        with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
            json.dump(dic_drugnames_DBID, f, ensure_ascii=False, indent=4)

    ######### Drug-Drug Interactions #########
    drug_drug_coodinates = []
    for drug_entry in tqdm(root):
        drug_drug_coodinates = get_drug_drug_coordinates(drug_entry, drug_drug_coodinates)
    #coordinates df
    df_d = pd.DataFrame(drug_drug_coodinates, columns=['Drug_ID_A', 'Drug_ID_B']) 
    df_d = df_d.drop_duplicates()
    logging.info(f'shape of coordinate file is: {df_d.shape}')
    logging.info('Writing coordinate file drug-drug effect...')
    df_d.to_csv(os.path.join(output_path, 'coordinates_drug_drug.tsv'), header=True,index=False ,sep="\t")





#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
