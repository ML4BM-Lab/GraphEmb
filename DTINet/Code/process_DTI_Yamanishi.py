import numpy as np
import pandas as pd
import random
import os, sys
from pubchempy import Compound
import requests
import logging
import xml.etree.ElementTree as ET
from tqdm import tqdm
import argparse
from re import search


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

#
def check_drug(drug_entry):
    try:
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        kegg_index = [ _.text for _ in drug_entry.findall('.//{http://www.drugbank.ca}resource')].index('KEGG Drug')
        kegg_ID = drug_entry.findall('.//{http://www.drugbank.ca}identifier')[kegg_index].text
        return drugbank_ID, kegg_ID
    except ValueError:
        return None, None

def get_dict_kegg2db():
    logging.debug('Reading DrugBank xml file...')
    tree = ET.parse('../../../Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    logging.debug('Succesfully read!')
    root = tree.getroot()
    dbids = []
    keggids = []
    for drug_entry in tqdm(root):
        drugbank_ID, kegg_ID = check_drug(drug_entry)
        if drugbank_ID and kegg_ID:
            dbids.append(drugbank_ID)
            keggids.append(kegg_ID)
    kegg2db = dict(zip(keggids, dbids))
    return kegg2db


def main():
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
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
        This script needs:
            - XXXXX
        
        Returns:
            - DTI tsv file
        '''
        )

    DB_PATH = args.dbPath
    logging.debug(f'DB_PATH=={DB_PATH}')
    logging.info(f'Working in output folder for: {DB_PATH}')
    db_name = get_DB_name(DB_PATH)
    check_and_create_folder(db_name)
    # Create relative output path
    wdir = os.path.join('../Data', db_name)
    logging.debug(f'working directory is: {wdir}')
    # wdir = '../Data/DrugBank'
    # load yamanishi data --> change when creating function <<<<< PARSE !
    # yamdb = 'E' # change here for go to other DB
    # wdir = '../Data/Yamanashi_et_al_GoldStandard/' + yamdb
    colnames_data = ['Kegg_ID', 'Gene']
    data_path = '../../../Data/Yamanashi_et_al_GoldStandard/' + DB_PATH.upper()+'/interactions/' + DB_PATH.lower() +'_admat_dgc_mat_2_line.txt'
    df = pd.read_csv(data_path, header = None, names = colnames_data, index_col=False, sep='\t')
    df['Gene'] = df['Gene'].map(lambda x: x[:3] + ":" + x[3:])
    logging.debug(f'{(df.head(2))}')

    # change from hsa: to Uniprot ---> wget http://rest.kegg.jp/conv/hsa/uniprot
    path_data_uniprot = '../Data/Yamanashi_et_al_GoldStandard/uniprot.txt'
    colnames_uniprot = ['Uniprot','Gene']
    df_uniprot = pd.read_csv(path_data_uniprot, header = None, names = colnames_uniprot, sep='\t')
    df_uniprot['Uniprot'] = df_uniprot['Uniprot'].map(lambda x: x.lstrip('up:'))
    hsa2uni = dict(zip(df_uniprot.Gene.tolist(), df_uniprot.Uniprot.tolist()))
    # repeated entries, but diff uniprots have the same sequence (so it's safe)
    df['Uniprot'] = df['Gene'].map(hsa2uni)
    logging.debug(f'{(df.head(2))}')
    # Change from KEGG Drug ID to DrugBank  
    logging.info("Loading kegg2db dict...")
    kegg2db = get_dict_kegg2db()
    df['DrugBank_ID'] = df['Kegg_ID'].map(kegg2db)
    # Process df before saving
    logging.info("Processing & saving DTIs to file...")
    df = df.drop(columns=['Gene', 'Kegg_ID']) # dejar al final cuando se caiga drug IDKegg
    df.columns = ['Protein', 'DrugBank_ID']
    DTI = df[['DrugBank_ID', 'Protein']] 
    # DTI.isna().sum()
    DTI = DTI.dropna()
    DTI.drop_duplicates()
    logging.debug(f'{(df.head(2))}')
    logging.info(f'Matrix Shape: {DTI.shape}')
    logging.info(f'unique drugs: {len(DTI.DrugBank_ID.unique())}')
    logging.info(f'unique targets: {len(DTI.Protein.unique())}')
    #DTI
    DTI.to_csv(os.path.join(wdir, f'DTI_{DB_PATH}.tsv'), header=True,index=False ,sep="\t")




#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################


