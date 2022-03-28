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

def get_compound_pubchem(drug):
	try:
		return Compound.from_cid(drug).isomeric_smiles
	except:
		return None

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
				#print("Hi" + smiles)
				#Chem.MolFromSmiles(str(smiles))
				break
	if not smiles:
		for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
			for ids in exids:
				if(ids.text == 'PubChem Substance'): 
					pubchem_id = exids[1].text
					smiles = get_compound_pubchem(pubchem_id)
					#print(smiles)
					break
	if not smiles:
		return(drugbank_ID, None)
	else:
		return(drugbank_ID, smiles)

tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/DrugBank/full_database.xml')
root = tree.getroot()
# list_drugs_w_smiles  = []
# for drug_entry in tqdm(root):
#list_drugs_w_smiles.append(get_list_drug_w_smiles(drug_entry))
list_drugs  = []
list_smiles = []
for drug_entry in tqdm(root):
	drug_id, smiles = get_smiles(drug_entry)
	if drug_id and smiles:
		list_drugs.append(drug_id)
		list_smiles.append(smiles) 
assert len(list_drugs) == len(list_smiles), 'The length of the Drug IDs does not match the number of SMILES'
dict_drugid_smiles = dict(zip(list_drugs, list_smiles))

df = pd.DataFrame(dict_drugid_smiles.items(), columns=['DrugBank_ID', 'SMILES'])

output_path = os.getcwd() + '/../Data/DrugBank/drugs_smiles.txt'
df.to_csv(output_path, header=None, index = None, sep = ' ')