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


'''
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
'''

logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)


wdir = '../Data/Davis_et_al'
db_file_path = '../../DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al.tsv'


davis = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence'])
davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
davis.head(4)



mod_davis = davis.copy()
mod_davis = mod_davis.drop(columns=['SMILES', 'GeneName', 'Sequence'])
##### work with drugs
mod_davis.loc[:, 'PubChemID'] = mod_davis.loc[:, 'PubChemID'].astype(int) # not float
mod_davis.loc[:, 'PubChemID'] = mod_davis.loc[:, 'PubChemID'].astype(str) # not float
mod_davis.head(2)

dic_cid_dbid = hf.pubchem_to_drugbankid()
mod_davis['DrugBankID'] = mod_davis['PubChemID'].map(dic_cid_dbid)
mod_davis.head()
# lens
drugs_before_request = len(mod_davis.DrugBankID.unique())
len(mod_davis.PubChemID.unique())

list_cids_wo_DB  = mod_davis[mod_davis['DrugBankID'].isna() == True].PubChemID.unique().tolist()

file_path_dict_cid2sid = os.path.join(wdir, 'dic_cid2sid.json')
cid2sid = hf.get_dic_cid2sid(list_cids_wo_DB)

# PUBCHEM SUBSTANCE
dic_kegg_db = hf.get_dict_kegg2db()
dic_kegg_pubchem = hf.get_dict_kegg2pubchemsid() # coge SUBSTANCE



#list_cids_wo_DB # CID
d1 ={} 
dic_pubchem2kegg = dict(zip(dic_kegg_pubchem.values(), dic_kegg_pubchem.keys()))
for drug_cid in list_cids_wo_DB:
    if drug_cid in cid2sid.keys():
        for element in cid2sid.get(drug_cid):
            if element in dic_kegg_pubchem.values():
                kegg = dic_pubchem2kegg.get(str(element))
                drugbank_id = dic_kegg_db.get(kegg)
                if drugbank_id:
                    d1[drug_cid] = drugbank_id

found_d1 = list(d1.keys())
logging.debug(f'found {len(found_d1)} drugs!')
list_cids_wo_DB = [x for x in list_cids_wo_DB if x not in found_d1]


# now with names
cid2syn = hf.get_dic_cid2syn(list_cids_wo_DB)

drugbank_dic = hf.drugname_drugbankid()
# make sure all is lower case for true comparison
drugbank_dic =  {k.lower(): v for k, v in drugbank_dic.items()}
cid2syn = {k: [i.lower() for i in v] for k, v in cid2syn.items()} 

d2 ={} 
for drug_cid in list_cids_wo_DB:
    if drug_cid in cid2syn.keys():
        for element in cid2syn.get(drug_cid):
            if element in drugbank_dic.keys():
                drugbank_id = drugbank_dic.get(element)
                if drugbank_id:
                    d2[drug_cid] = drugbank_id


### jj
mod_davis['DrugBankID_A'] = mod_davis['PubChemID'].map(d1)
logging.debug( len(mod_davis['DrugBankID'].unique()) )

mod_davis['DrugBankID_B'] = mod_davis['PubChemID'].map(d2)
logging.debug(len(mod_davis['DrugBankID_B'].unique()))

mod_davis.DrugBankID.fillna(mod_davis.DrugBankID_A, inplace=True)
mod_davis.DrugBankID.fillna(mod_davis.DrugBankID_B, inplace=True)

drugs_after_resquest = len(mod_davis.DrugBankID.unique())
logging.debug(f'we did our best and rescued {drugs_after_resquest-drugs_before_request} drugs')




