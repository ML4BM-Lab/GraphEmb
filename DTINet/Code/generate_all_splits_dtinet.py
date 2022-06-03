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
    script for generating all splits for a single folder 
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
    parser.add_argument("dbPath", help="Select a Dataset!", type=str)
    args = parser.parse_args()
    DB_PATH = args.dbPath
    script_spits = 'generate_splits_dtinet.py'
    splits = ['Sp', 'Sd', 'St']
    subsamplings = [True, False]
    for sub in subsamplings:
        for split_type in splits:
            try:
                if sub:
                    logging.info(f'Generating {split_type} splits for {DB_PATH} with subsampling')
                    return_code = sp.check_call(['python3', script_spits, '--dbPath', DB_PATH, '--split_type', split_type, '-subsampling'])
                else: 
                    logging.info(f'Generating {split_type} splits for {DB_PATH} without subsampling')
                    return_code = sp.check_call(['python3', script_spits, '--dbPath', DB_PATH, '--split_type', split_type])
                # check return code
                if return_code ==0: 
                    logging.info('EXIT CODE 0')
            except sp.CalledProcessError as e:
                logging.info(e.output)
                break

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
