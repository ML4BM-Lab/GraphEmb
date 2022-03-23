import os, sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import multiprocessing as mp
import time
import json
from pubchempy import Compound
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
import xml.etree.ElementTree as ET
from itertools import repeat
import requests
import logging
from re import search
import argparse
import argparse


log_levels = {
    0: logging.CRITICAL,
    1: logging.ERROR,
    2: logging.WARN,
    3: logging.INFO,
    4: logging.DEBUG,
}
logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)

###############
def get_compound_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles


def get_drug_nodes(file_path_drugs):
	# need to define al global the drug matrix !!
	global drug_drug, drug_dis, drug_se, DTI
	# check
	if (not os.path.exists(file_path_drugs)):
		logging.info('Getting Drug Nodes...')
		## get a list of tuples (drug, smiles)

		logging.info('Reading DrugBank xml file')
		tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/DrugBank/full_database.xml')
		root = tree.getroot()
		logging.info('Retrieving a list of drugs with SMILES')
		list_drugs_w_smiles  = []
		for drug_entry in tqdm(root):
			#list_drugs_w_smiles.append(get_list_drug_w_smiles(drug_entry))
		# 
		drugs_linked_to_drugs = set(drug_drug.D1.tolist() + drug_drug.D2.tolist()) # 4418
		drugs_linked_to_disease = set(drug_dis.DrugBankID.tolist()) # 2754
		drugs_linked_to_sideeffect = set(drug_se.DrugBank_ID.tolist()) # 675
		drugs_linked_to_proteins = set(DTI.Drug.tolist()) # 7627 ###### <--------- THIS IS THE ONLY ONE THAT CHANGES FOR EACH DATABASE
		not_isolated_drugs = list(drugs_linked_to_drugs.union(drugs_linked_to_disease, drugs_linked_to_sideeffect, drugs_linked_to_proteins )) ### 9720
		not_isolated_drugs.sort()
		# intersec with available SMILES
		list_of_drug_nodes = list(set(not_isolated_drugs).intersection(set(list_drugs_w_smiles))) #  ##############---->>> ** aqui ya estaria!
		list_of_drug_nodes.sort()
		#len(set(not_isolated_drugs).intersection(set(list_drug_w_smiles))) # 8638
		logging.info(f'There are {len(list_of_drug_nodes)} isolated drugs with available smiles from {len(not_isolated_drugs)} total isolated nodes (diff: {len(not_isolated_drugs)-len(list_of_drug_nodes)})')
		np.savetxt(os.path.join(file_path_drugs), list_of_drug_nodes, fmt='%s')
	else:
		list_of_drug_nodes = np.loadtxt(file_path_drugs, dtype='str').tolist()
	#
	return list_of_drug_nodes


def get_smiles(drug_entry):
	'''
	Get list of drugs with smiles from DrugBank xml, 
	for SMILES that are not in DrugBank retrieves them from PubChem
	Then, check if we can create a fingerprint (if not, do not include)
	'''
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	smiles = None
	fp = None
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
				break
	if not smiles:
		for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
			for ids in exids:
				if(ids.text == 'PubChem Compound'): 
					pubchem_id = exids[1].text
					smiles = get_compound_pubchem(pubchem_id)
					break
	if not smiles:
		return(drugbank_ID, None)
	elif Chem.MolFromSmiles(str(smiles)):
		return(drugbank_ID, smiles)
	return(drugbank_ID, None)

	""" if smiles:
		# check if fp exists
		fp = Chem.MolFromSmiles(str(smiles))
		if fp:
			#drug = (drugbank_ID)
			return (drugbank_ID, smiles)
		else: 
			return (drugbank_ID, None)
	else: 
		return (drugbank_ID, smiles)

 """


logging.info('Reading DrugBank xml file')
tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/DrugBank/full_database.xml')

root = tree.getroot()

list_drugs  = []
list_smiles = []
for drug_entry in tqdm(root):
	drug_id, smiles = get_smiles(drug_entry)
	if drug_id and smiles:
		list_drugs.append(drug_id)
		list_smiles.append(smiles)

assert len(list_drugs) == len(list_smiles), 'The length of the Drug IDs does not match the number of SMILES'

	