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
import helper_functions_dtinet as hf




######################################## START MAIN #########################################
#############################################################################################

def main():
	'''
	Seq run:
	run get_coord for BIOSNAP
	DTI already as in DB folder
	run get_all matrix BIOSNAP

	'''
	parser = argparse.ArgumentParser() 
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
					help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	parser.add_argument("dbPath", help="Path to the database output 'DrugBank' ", type=str)
	args = parser.parse_args()
	DB_PATH = args.dbPath
	if DB_PATH != 'BIOSNAP':
		raise NameError('This script is for BIOSNAP  Dabatase')
	# check folder
	db_name = hf.get_DB_name(DB_PATH)
	hf.check_and_create_folder(db_name)
	# list of scripts to execute sequentially
	list_of_pys = ['get_coord.py', 'get_all_matrix_BIOSNAP.py']
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
	# esto quitarlo de aqui y que get coord sea general para todos 

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################

