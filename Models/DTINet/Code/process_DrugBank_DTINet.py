import os, sys
import argparse
import logging
import pandas as pd
import json
from tqdm import tqdm
from re import search
import xml.etree.ElementTree as ET
import helper_functions_dtinet as hf




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
    logging.info("============== DrugBank ==============")
    ###  output details
    logging.info(
        '''
        This script needs:
            - DrugBank xml file

        This script generates the following files:
            - all_drugs.txt
            - all_drug_dic_map.txt (/ dic_drugnames_DBID.json)
            - edgelist_drug_drug.tsv
        
        DTINet uses from DrugBank the processed files:
            - mat_drug_drug.txt
        '''
        )

    # OUTPUT DIRECTORY
    DB_PATH = args.dbPath
    db_name = hf.get_DB_name(DB_PATH)
    hf.check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)

    # info 
    logging.info(f'Processing DrugBank in {output_path}')

    #
    logging.debug('Reading DrugBank xml file...')
    tree = ET.parse('../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    root = tree.getroot()

    ########### DRUG NODES ###########
    drug_IDs = []
    drug_names = []
    for drug_entry in tqdm(root, desc='Retrieving Drug IDs from DrugBank', position=0, leave=True):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        name = drug_entry.find('{http://www.drugbank.ca}name').text 
        drug_names.append(name)
        drug_IDs.append(drugbank_ID)

    # Drug nodes (all_drugs.txt). All nodes in DB (w\o filter) 
    logging.debug(f'Writing all_drugs.txt...')
    with open(os.path.join(output_path, 'all_drugs.txt'), 'w') as f:
        for item in drug_IDs:
            _ = f.write("%s\n" % item)
    
    # Write dictionary (drugbank_ID and name)
    logging.debug(f'Writing all_drug_dic_map.txt...')
    with open(os.path.join(output_path,'all_drug_dic_map.txt'), 'w') as f:
        for i in range(len(drug_IDs)):
            _ = f.write("%s:%s\n" % (drug_IDs[i],drug_names[i]))
    # json
    if args.json == True:
        dic_drugnames_DBID = dict(zip(drug_IDs, drug_names)) 
        file_name_json = 'dic_drugnames_DBID.json'
        logging.debug(f'Writing {file_name_json}...')
        with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
            json.dump(dic_drugnames_DBID, f, ensure_ascii=False, indent=4)

    ######### Drug-Drug Interactions #########
    logging.info(f'    Getting Drug-Drug coordinates...')
    drug_drug_coodinates = []
    for drug_entry in tqdm(root, desc='Retrieving drug-drug interactions from DrugBank',  position=0, leave=True):
        drug_drug_coodinates = hf.get_drug_drug_coordinates(drug_entry, drug_drug_coodinates)
    
    #coordinates df
    df_d = pd.DataFrame(drug_drug_coodinates, columns=['Drug_ID_A', 'Drug_ID_B']) 
    df_d = df_d.drop_duplicates()
    logging.info(f'    shape of coordinate file is: {df_d.shape}')
    logging.debug('Writing coordinate file drug-drug effect...')
    df_d.to_csv(os.path.join(output_path, 'edgelist_drug_drug.tsv'), header=True,index=False ,sep="\t")





#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####------------------------------------------------------------------------------------------
####################### END OF THE CODE ########################################################
