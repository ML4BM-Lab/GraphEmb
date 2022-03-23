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
#### 

'''
This script should execute:
    - process_HPRD_DTINet.py
    - process_SIDER_DTINet.py
    - process_CTD_DTINet.py
    - process_DrugBank.py

=> parse argument for DB directory 
'''


######################################## START MAIN #########################################
#############################################################################################

def main():
    '''
    This executes 4 complementary scripts
    and writes 5 coordinate files:
        - coordinates_PPI.tsv
        - coordinates_drug_se.tsv
        - coordinates_protein_disease.tsv
        - coordinates_drug_disease.tsv
        - coordinates_drug_drug.tsv

    The only left is DTI, this changes for each model changes for each model
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
    parser.add_argument("dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)

    args = parser.parse_args()

    DB_PATH = args.dbPath

    list_of_pys = ['process_HPRD_DTINet.py', 'process_SIDER_DTINet.py', 'process_CTD_DTINet.py', 'process_DrugBank.py']
    
    # check that exception works; check other options
    for script in list_of_pys:
        try:
            return_code = sp.check_call(['python3', script, DB_PATH])
            if return_code ==0: 
                print('EXIT CODE 0')
        except sp.CalledProcessError as e:
            print(e.output)
            break

    if DB_PATH in ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank', 'E', 'GPCR', 'IC', 'NR'):
        try:
            if DB_PATH == 'DrugBank':
                print('Already as tsv')
            else:
                return_code = sp.check_call(['python3', f'process_DTI_{DB_PATH}_new.py', DB_PATH])
            #
            if return_code == 0:
                print('EXIT CODE 0')
        except sp.CalledProcessError as e:
            print(e.output)
        except FileNotFoundError:
            print(f"No 'process_DTI_{DB_PATH}_new.py' file yet")
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
