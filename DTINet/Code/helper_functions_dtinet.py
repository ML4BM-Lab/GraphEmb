import os, sys, uuid
import argparse
import logging
import numpy as np
import pandas as pd
import json
import requests
from tqdm import tqdm
from re import search
import time
import xml.etree.ElementTree as ET
import pubchempy as pcp
from pubchempy import Compound
import subprocess as sp
from pubchempy import Compound
from rdkit import Chem
from rdkit import DataStructs

'''
from importlib import reload
reload(hf)
'''

# Databases Paths
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

# funct
def get_cid(compound_name):
    try:
        comp_id = pcp.get_compounds(str(compound_name), 'name')[0].cid
        time.sleep(1) # increase if error
    except IndexError:
        comp_id = np.nan
    return comp_id

def get_compound_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

def pubchem_to_drugbankid():
    '''
    Parse  DrugBank DB xml and 
    create a dictionary from PubChemID to DrugBankID
    '''
    logging.debug('Reading DrugBank xml database')
    tree = ET.parse('../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    logging.debug('Succesfully read!')
    root = tree.getroot()
    dbids = []
    pubchemids = []
    for drug_entry in tqdm(root, desc='Retrieving Drugbank  & PubChem ID', position=0, leave=True):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        pubchem_id = np.nan # to not repeat id if it does not appear later
        for props in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
            for prop in props:
                if(prop.text == 'PubChem Compound'): 
                    pubchem_id = props[1].text
                break # romper una vez que encuentre el pubchem 
        dbids.append(drugbank_ID)
        pubchemids.append(pubchem_id)
    #
    dic_cid_dbid = dict((zip(pubchemids, dbids)))
    # first element is nan, delete
    dic_cid_dbid[list(dic_cid_dbid.keys())[0]]
    dic_cid_dbid.pop(list(dic_cid_dbid.keys())[0])
    return dic_cid_dbid

def get_drug_drug_coordinates(drug_entry, drug_drug_coodinates):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	all_interactions =  drug_entry.findall('.//{http://www.drugbank.ca}drug-interactions')[0]
	for inter in all_interactions: # aqui iteramos en targets
		coordinate = (drugbank_ID, list(inter)[0].text)	
		drug_drug_coodinates.append(coordinate)
	return drug_drug_coodinates

def select_drugbankid(x):
    '''
    This function returns the pubchem id depending on what code is available in data
	in order to retrieve the max possible information.
    '''
    if x['Pubchem_map_flat'] == x['Pubchem_map_stereo']:
        return x['Pubchem_map_flat']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and (x['Pubchem_map_stereo']) != np.nan ) :
        return x['Pubchem_map_stereo']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and (np.isnan(x['Pubchem_map_flat']) == False)):
        return x['Pubchem_map_flat']
    else:
        return np.nan

def read_fasta(path):
	names=[]
	seqs = []
	with open(path, 'r') as f:
		for line in f:
			if line.startswith('>'):
				names.append(line.strip().replace('>', ''))
			else:
				seqs.append(line.strip())
	return zip(names, seqs)

def write_fasta(path, target, seq):
	fl_name = os.path.join(path, target.replace(':', '_')+'.fasta')
	if os.path.exists(fl_name):
		logging.debug(f'File {fl_name} already exists')
		return fl_name
	with open(fl_name, 'w') as f:
		_ = f.write('>'+target+'\n'+seq+'\n')
	return fl_name

def create_remove_tmp_folder(path):
	if not os.path.exists(path):
		logging.info('Creating tmp folder: {}'.format(path))
		os.makedirs(path)
		return path
	else: 
		return path


def check_and_create_fasta(target, seq):
	global PATH
	fasta1 = os.path.join(PATH, target.replace(':', '_')+'.fasta')
	if not os.path.exists(fasta1):
		fasta1 = write_fasta(PATH, target, seq)
	return fasta1

# SW

def get_amino_uniprot(proteinID):
    r = requests.get(f'https://www.uniprot.org/uniprot/{proteinID}.fasta')
    if r.status_code == 200 and r.text:
        return (proteinID, ''.join(r.text.split('\n')[1:]))
    else: 
        #print('Protein sequence not found in uniprot database')
        return (proteinID, None)


def get_pairwise_tanimoto(smiles1,smiles2, dic): #dic == smile2fp
	try:
		for smile in [smiles1, smiles2]:
			if not smile in dic:
				mol1 = Chem.MolFromSmiles(str(smile))
				fp1  = Chem.RDKFingerprint(mol1)
				dic[smile] = fp1
		tani = DataStructs.FingerprintSimilarity(dic[smiles1],dic[smiles2]) #pairwise similarity
		return tani
	except:
		return None


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
    tree = ET.parse('../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
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


def get_dict_hsa2uni():
	'''
	http://rest.kegg.jp/conv/hsa/uniprot 
	return a dict to change from hsa to uniprot ID
	'''
	r = requests.get('http://rest.kegg.jp/conv/hsa/uniprot')
	rlines = r.text.split('\n')
	up, hsa = [], []
	for line in rlines:
		if line != '':
			upi, hsai = line.split('\t')
			up.append(upi.lstrip('up:'))
			hsa.append(hsai)
	return dict(zip(hsa, up))

# new from get_SW_score_Yamanishi.py in DTI2Vec
def get_SW_score(pair1, pair2, tmp_path):
	target1, _ = pair1
	target2, _ = pair2
	fasta1 = os.path.join(tmp_path, target1.replace(':', '_')+'.fasta')
	fasta2 = os.path.join(tmp_path, target2.replace(':', '_')+'.fasta')
	args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
			'-asequence', fasta1 , '-bsequence', fasta2, 
			'-gapopen', '10.0', '-gapext', '0.5', 
			'-stdout']
	try:
		score = sp.Popen(args, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.DEVNULL)
		score = score.communicate(b"\n")
		score = score[0].decode().split('\n')
		score = extract_score(score)
		return score
	except:
		logging.warning(f'Not able to compute SW score for : {target1}, {target2}')

def extract_score(score):
	score = [line for line in score if line.startswith('# Score:')]
	if score:
		return float(score[0].split()[-1])
	else:
		return None

def write_all_fastas(fastas, path):
	for header, seq in fastas:
		write_fasta(path, header, seq)
	logging.info('All fastas written')
