from ast import keyword
import requests
from requests.adapters import HTTPAdapter, Retry
import re
import logging
import numpy as np
import pandas as pd
import time
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

###

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)
###


def get_key(ecode, kws):
    global class_list
    class_list = np.loadtxt('../Data/keys.txt', dtype='str', delimiter='\n').tolist()
    class_list = [key.lower() for key in class_list]
    logging.debug("""keywords from: https://www.uniprot.org/keywords/KW-9992; 
                    deleting enzyme keywords because using EC Code for those""")
    if ecode:
        mol_funct = 'enzyme'
    elif any(key in class_list for key in kws):
        mol_funct = [word for word in class_list if word in kws][0]
    else:
        mol_funct = 'notAnnotated' 
    return mol_funct


def get_json_info_uni(data):
    check_uni_id = data['primaryAccession']
    length = int(data['sequence']['length'])
    ecode = None
    if data.get('proteinDescription').get('recommendedName'):
        recname = data.get('proteinDescription').get('recommendedName')
        if recname.get('ecNumbers'):
            ecode = recname.get('ecNumbers')[0].get('value')
    mol_funct = None
    if data.get('keywords'):
        dic_molfunc = [line for line in data.get('keywords') if line.get('category') == 'Molecular function']
        keywords = [line.get('name').lower() for line in dic_molfunc]
        mol_funct = get_key(ecode, keywords) 
    res = (check_uni_id, length, ecode, mol_funct)
    return res



# list of proteins w structure
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins_struct = list_pdbs + list_afold


PKLS_PATH = '../Data/pkls'
prot_rmsd = pd.read_pickle(os.path.join(PKLS_PATH, 'RMSD_full_matrix.pkl')) #
list_proteins = prot_rmsd.index.tolist()

#list_proteins =  list_proteins_struct # np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()
n = 500
list_proteins = [list_proteins[i:i + n] for i in range(0, len(list_proteins), n)]

full_info= []
for lts in tqdm(list_proteins, desc='Retrieving uniprot annotation'):
    unilist = ','.join(lts)
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/accessions?accessions={unilist}')
    jsons = r.json()['results']
    for data in jsons:
        res = get_json_info_uni(data)
        full_info.append(res)


cols_all = ['UniprotID', 'length', 'EC', 'mol_func']
df_all = pd.DataFrame.from_records(full_info, columns=cols_all)

lost_proteins = list(set(list_proteins_struct) - set(df_all.UniprotID))
logging.debug(f'{len(lost_proteins)} proteins changed un uniprot, accesing secondary')

# those can be queried with the secondary accesion number 
for prot in tqdm(lost_proteins, desc='Secondary Acc'):
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/search?&query=sec_acc:{prot}')
    data = r.json()['results'][0]
    res=get_json_info_uni(data)
    df_all.loc[len(df_all)] =  res



# clear enzyme
df_all['EC'][~df_all.EC.isna()] = df_all.EC[~df_all.EC.isna()].apply(lambda x: x[:3])

# Add Source
dic_alfa = dict(zip(list_afold, ['AFold'] * len(list_afold)))
dict_pdb = dict(zip(list_pdbs, ['PDB'] * len(list_pdbs)))
dic_source = {}
dic_source.update(dic_alfa)
dic_source.update(dict_pdb)
df_all['source'] = df_all.UniprotID.map(dic_source)

df_all[df_all.source.isna()]


# changes
# {'P06312': 'P01911',# 'Q8ILQ7': 'P01911'}

df_all  = df_all[~df_all.source.isna()]

OUT_FILE = '../Data/pkls/df_annot_proteins_structure.pkl'
df_all.to_pickle(OUT_FILE)

logging.debug(df_all.head())

#exit(0)

set(prot_rmsd.index.tolist()) - set(df_all.UniprotID.tolist()) 
len(prot_rmsd.index.unique().tolist())
len(df_all.UniprotID.unique().tolist())


Annot_Prot = pd.DataFrame(prot_rmsd.index.unique().tolist(), columns=['UniprotID'])

final_annot = pd.merge(Annot_Prot, df_all, on='UniprotID', how='right').fillna('-')

final_annot.to_pickle('../Data/pkls/final_protein_annot.pkl')