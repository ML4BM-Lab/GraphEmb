
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


def select_drugbankid(x):
    '''
    *********************** explain  *******
    '''
    if x['Pubchem_map_flat'] == x['Pubchem_map_stereo']:
        return x['Pubchem_map_flat']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and ( x['Pubchem_map_stereo']) != np.nan ) :
        return x['Pubchem_map_stereo']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and (np.isnan(x['Pubchem_map_flat']) == False)):
        return x['Pubchem_map_flat']
    else:
        return np.nan

def drugbankid_to_pubchem():
    logging.debug('Reading database')
    tree = ET.parse('../../../Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    logging.debug('Succesfully read!')
    root = tree.getroot()
    dbids = []
    pubchemids = []
    for drug_entry in tqdm(root):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        pubchem_id = np.nan # to not repeat id if it does not appear later
        for props in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
            for prop in props:
                if(prop.text == 'PubChem Compound'): 
                    pubchem_id = props[1].text
                break # romper una vez que encuentre el pubchem 
        dbids.append(drugbank_ID)
        pubchemids.append(pubchem_id)
    #
    dic_cid_dbid = dict((zip(pubchemids, dbids)))
    # first element is nan, delete
    dic_cid_dbid[list(dic_cid_dbid.keys())[0]]
    dic_cid_dbid.pop(list(dic_cid_dbid.keys())[0])
    return dic_cid_dbid

######################################## START MAIN #########################################
#############################################################################################

def main():
    '''
    SIDER 
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
    #
    # OUTPUT DIRECTORY
    DB_PATH = args.dbPath
    db_name = get_DB_name(DB_PATH)
    check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)
    ## SIDER DATA FOLDER
    data_path = '../../../Data/cross_side_information_DB/SIDER'
    # info 
    logging.info(f'Processing {data_path[40:]} in {output_path}')
    #
    ####### read DrugBank
    logging.info('Reading DrugBank xml file...')
    dic_cid_dbid = drugbankid_to_pubchem()
    #
    ####### DRUG LIST
    file_all_se = 'meddra_all_se.tsv'
    logging.info(f'Using: {file_all_se}.csv')
    df_drug_se = pd.read_csv(os.path.join(data_path, file_all_se), sep='\t',
                names = ['STITCH_flat', 'STITCH_stereo', 'x1', 'x2', 'x3', 'se'],
                usecols = ['se', 'STITCH_flat', 'STITCH_stereo'])
    # lowercase
    df_drug_se.se = df_drug_se.iloc[:, 2:3].applymap(lambda s: s.lower() if type(s) == str else s)
    se_list = df_drug_se.se.unique().tolist()
    # save se.txt
    logging.info(f'Writing se.txt')
    np.savetxt(os.path.join(output_path ,'se.txt'), se_list , newline='\n', fmt='%s')
    #
    ####### COORDINATES DRUG - SE: Preproc
    logging.info('Processing data...')
    ## take both to get the max posible number of DBIDs
    df_drug_se['Pubchem_flat'] = df_drug_se['STITCH_flat'].apply(lambda x: x[4:]).astype(int)
    df_drug_se['Pubchem_stereo'] = df_drug_se['STITCH_stereo'].apply(lambda x: x[4:]).astype(int)
    # make sure string 
    df_drug_se['Pubchem_flat'] = df_drug_se['Pubchem_flat'].astype(str)
    df_drug_se['Pubchem_stereo'] = df_drug_se['Pubchem_stereo'].astype(str)
    # map from pubchem to DrugBank ID
    df_drug_se['Pubchem_map_flat'] = df_drug_se['Pubchem_flat'].map(dic_cid_dbid)
    df_drug_se['Pubchem_map_stereo'] = df_drug_se['Pubchem_stereo'].map(dic_cid_dbid)
    # We only need to columns, for safety deep copy 
    data_drug_to_se = df_drug_se[['se', 'Pubchem_map_flat', 'Pubchem_map_stereo']].copy(deep=True)
    # Create the column DrugBank_ID applying the function 
    data_drug_to_se['DrugBank_ID'] = data_drug_to_se.apply(select_drugbankid, axis=1)
    #  not take only the interesting columns & drop NaN & duplicates
    coordinates_drug_se = data_drug_to_se[['DrugBank_ID', 'se']]
    coordinates_drug_se = coordinates_drug_se.dropna()
    coordinates_drug_se = coordinates_drug_se.drop_duplicates()
    # finish
    logging.info(f'shape of coordinate file is: {coordinates_drug_se.shape}')
    logging.info('Writing coordinate file drug-side effect...')
    coordinates_drug_se.to_csv(os.path.join(output_path, 'coordinates_drug_se.tsv'), header=True, index=False ,sep="\t")



#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
