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
import subprocess as sp
from rdkit import RDLogger                     



######################################## START MAIN #########################################
#############################################################################################

def main():
    '''
    This executes 4 complementary scripts
    and writes 5 coordinate files:
        - edgelist_PPI.tsv
        - edgelist_drug_se.tsv
        - edgelist_protein_disease.tsv
        - edgelist_drug_disease.tsv
        - edgelist_drug_drug.tsv

    The only left is DTI, this changes for each model changes for each model
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
    parser.add_argument("dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')  # disable RDKit Log
    DB_PATH = args.dbPath
    list_of_pys = ['process_HPRD_DTINet.py', 'process_SIDER_DTINet.py', 'process_CTD_DTINet.py', 'process_DrugBank_DTINet.py']
    # check that exception works; check other options
    for script in list_of_pys:
        try:
            #logging.info(f'Running {script}')
            return_code = sp.check_call(['python3', script, DB_PATH])
            if return_code ==0: 
                logging.info('EXIT CODE 0')
        except sp.CalledProcessError as e:
            logging.info(e.output)
            break

    logging.info(
        '''
        ------------------- FINISHED: get_coord.py ---------------------
        '''
        )

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################

