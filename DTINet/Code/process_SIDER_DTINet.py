import os, sys
import argparse
import logging
import numpy as np
import pandas as pd
import json
import requests
from re import search
import xml.etree.ElementTree as ET
import helper_functions_dtinet as hf




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
    logging.info("============== SIDER ==============")
    # Output details
    logging.info(
        '''
        This script needs:
            - From SIDER: meddra_all_se.tsv, 
            - DrugBank xml file
        This script generates the following files:
            - se.txt
            - edgelist_drug_se.tsv
        DTINet uses from SIDER the processed files:
            - mat_drug_disease.txt
            - mat_protein_disease.txt
        '''
        )
    #
    # OUTPUT DIRECTORY
    DB_PATH = args.dbPath
    db_name = hf.get_DB_name(DB_PATH)
    hf.check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)
    ## SIDER DATA FOLDER
    data_path = '../../DB/Data/cross_side_information_DB/SIDER'
    logging.info(f'Processing {data_path[40:]} in {output_path}')
    ## Read DrugBank
    logging.debug('Reading DrugBank xml file...')
    dic_cid_dbid = hf.pubchem_to_drugbankid()

    ## DRUG LIST
    file_all_se = 'meddra_all_se.tsv'
    logging.debug(f'Using: {file_all_se}.csv')
    df_drug_se = pd.read_csv(os.path.join(data_path, file_all_se), sep='\t',
                names = ['STITCH_flat', 'STITCH_stereo', 'x1', 'x2', 'x3', 'se'],
                usecols = ['se', 'STITCH_flat', 'STITCH_stereo'])
    # lowercase
    df_drug_se.se = df_drug_se.iloc[:, 2:3].applymap(lambda s: s.lower() if type(s) == str else s)
    se_list = df_drug_se.se.unique().tolist()
    # save se.txt
    logging.debug(f'Writing se.txt...')
    np.savetxt(os.path.join(output_path ,'se.txt'), se_list , newline='\n', fmt='%s')

    ####### COORDINATES DRUG - SE: Preproc
    logging.info(f'    Getting Drug-Se coordinates...')
    logging.debug('Processing data...')
    ## take both to get the max posible number of DBIDs
    #df_drug_se['Pubchem_flat'] = df_drug_se['STITCH_flat'].apply(lambda x: x[4:]).astype(int)
    #df_drug_se['Pubchem_stereo'] = df_drug_se['STITCH_stereo'].apply(lambda x: x[4:]).astype(int)
    regex = r'CID[\d]{1}'
    df_drug_se['Pubchem_flat'] = df_drug_se.STITCH_flat.replace(regex, '', regex=True).astype(int)
    df_drug_se['Pubchem_stereo'] = df_drug_se.STITCH_stereo.replace(regex, '', regex=True).astype(int)
    # make sure string 
    df_drug_se['Pubchem_flat'] = df_drug_se['Pubchem_flat'].astype(str)
    df_drug_se['Pubchem_stereo'] = df_drug_se['Pubchem_stereo'].astype(str)
    # map from pubchem to DrugBank ID
    df_drug_se['Pubchem_map_flat'] = df_drug_se['Pubchem_flat'].map(dic_cid_dbid)
    df_drug_se['Pubchem_map_stereo'] = df_drug_se['Pubchem_stereo'].map(dic_cid_dbid)
    # We only need to columns, for safety deep copy 
    data_drug_to_se = df_drug_se[['se', 'Pubchem_map_flat', 'Pubchem_map_stereo']].copy(deep=True)
    # Create the column DrugBank_ID applying the function 
    data_drug_to_se['DrugBank_ID'] = data_drug_to_se.apply(hf.select_drugbankid, axis=1)
    #  not take only the interesting columns & drop NaN & duplicates
    coordinates_drug_se = data_drug_to_se[['DrugBank_ID', 'se']]
    coordinates_drug_se = coordinates_drug_se.dropna()
    coordinates_drug_se = coordinates_drug_se.drop_duplicates()
    # finish
    logging.info(f'    Shape of coordinate file is: {coordinates_drug_se.shape}')
    logging.debug('Writing coordinate file drug-side effect...')
    coordinates_drug_se.to_csv(os.path.join(output_path, 'edgelist_drug_se.tsv'), header=True, index=False ,sep="\t")


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####------------------------------------------------------------------------------------------
####################### END OF THE CODE #######################################################
