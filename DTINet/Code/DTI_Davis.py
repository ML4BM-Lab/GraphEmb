import os
import numpy as np
import pandas as pd
import logging
import argparse
import helper_functions_dtinet as hf




######################################## START MAIN #########################################
#############################################################################################

def main():
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=4,
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
    logging.info("============== DTIs Yamanishi ==============")
    ### log output detals
    logging.info(
        '''
        This script needs:
            - Drug Target Interactions for Davis et al
        Returns:
            - DTI tsv file --- only? 
        '''
        )

    DB_PATH = args.dbPath
    logging.debug(f'DB_PATH=={DB_PATH}')
    logging.info(f'Working in output folder for: {DB_PATH}')
    db_name = hf.get_DB_name(DB_PATH)
    hf.check_and_create_folder(db_name)
    # Create relative output path
    wdir = os.path.join('../Data', db_name)
    logging.debug(f'working directory is: {wdir}')
    

db_file_path = '../../DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al.tsv'
davis = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence'])
davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
davis.head(4)

logging.debug("get all needed dictionaries")

# Pubchem to drugbank id parsing DrugBank
dic_cid_dbid = hf.pubchem_to_drugbankid()



mod_davis = davis.copy()
# mod_davis = mod_davis.drop(columns=['SMILES', 'GeneName', 'Sequence'])

##### work with drugs
mod_davis.loc[:, 'PubChemID'] = mod_davis.loc[:, 'PubChemID'].astype(int) # not float
mod_davis.loc[:, 'PubChemID'] = mod_davis.loc[:, 'PubChemID'].astype(str) # not float
mod_davis.head(2)


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-----------------------------------------------------------------------------------------
####################### END OF THE CODE ######################################################

