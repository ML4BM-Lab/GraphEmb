import pandas as pd
from tqdm import tqdm
import helper_functions_dtinet as hf
import pubchempy as pcp
from tqdm import tqdm
import os
import json
import time
import random
import requests
import logging


logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)


def get_cid2sid(file_path_dict_cid2sid, list_cids):
	if (not os.path.exists(file_path_dict_cid2sid)):
		logging.debug('dic_cid2sid.json does not exist, creating it...')
		cid2sid = hf.get_dic_cid2sid(list_cids)
		with open(file_path_dict_cid2sid, 'w') as outfile:
			json.dump(cid2sid, outfile)
	else:
		logging.debug('Reading dic_cid2sid.json')
		with open(file_path_dict_cid2sid, 'r') as f:
			cid2sid = json.load(f)
	return cid2sid




wdir = '../Data/BindingDB'
db_file_path = '../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
# load 
bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence'])
bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID', 'Target Sequence': 'Sequence'}, axis=1)

#### for proteins
list_of_protein_nodes = bindingdb.UniprotID.unique().tolist()
# estas luego se preprocesan como lo otro
dict_protein_sequence = dict(zip(bindingdb.UniprotID, bindingdb['Sequence'])) # done just name as in all_mat

####### for drugs ########## ---> working here! 
mod_bindingdb = bindingdb.copy()
mod_bindingdb = mod_bindingdb.drop(columns=['SMILES', 'UniprotID', 'Sequence'])
##### work with drugs
mod_bindingdb.loc[:, 'PubChemID'] = mod_bindingdb.loc[:, 'PubChemID'].astype(int) # not float
mod_bindingdb.loc[:, 'PubChemID'] = mod_bindingdb.loc[:, 'PubChemID'].astype(str) # not float
mod_bindingdb.head(2)

dic_cid_dbid = hf.pubchem_to_drugbankid()
mod_bindingdb['DrugBankID'] = mod_bindingdb['PubChemID'].map(dic_cid_dbid)
mod_bindingdb.head(6)

list_cids_wo_DB  = mod_bindingdb[mod_bindingdb['DrugBankID'].isna() == True].PubChemID.unique().tolist()
len(list_cids_wo_DB)




file_path_dict_cid2sid = os.path.join(wdir, 'dic_cid2sid.json')

cid2sid = get_cid2sid(file_path_dict_cid2sid, list_cids_wo_DB)
 

'''
# other read file from json

#request to 

# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/51/sids/TXT

# dict to create 
# cid2sid 

# from pubchem to drugbank
dic_kegg_db = hf.get_dict_kegg2db()
dic_kegg_pubchem = hf.get_dict_kegg2pubchem() # coge SUBSTrANCE

dic_pubchem2kegg = dict(zip(dic_kegg_pubchem.values(), dic_kegg_pubchem.keys()))


bindingdb['KEGG'] = bindingdb['PubChemID'].map(dic_pubchem2kegg)


dtis = bindingdb[['DrugBankID', 'UniprotID', 'SMILES', 'Target Sequence']]
dtis.shape
dtis = dtis.dropna() # 25699 dtis lost

dtis.shape

len(dtis.DrugBankID.unique())
len(dtis.UniprotID.unique())

#dtis = dtis.dropna() # 1513dro
#dtis[dtis.Label == 1]


##### create a dict (for later)
from rdkit import Chem


list_of_drug_nodes = bindingdb.DrugBankID.unique().tolist()
# check if smile is ok
# maybe better as apply fp in dataset
def check_fp(smiles):
	a = Chem.MolFromSmiles(str(smiles))
	# probar con fp lo mismo? 
	return a

dtis['SMILES'].apply(lambda x: Chem.MolFromSmiles(str(x)))

dict_drugid_smiles = dict(zip(bindingdb.DrugBankID, bindingdb.SMILES)) 
# but need to check if tanimoto can be applied !  => fp!!!


'''