import os
from traceback import print_tb
import numpy as np
import pandas as pd
import glob
import logging
import requests
from tqdm import tqdm
import pandas as pd
import os
import logging
import re
import csv


# funct 
def get_key(df_line, EC, class_list):
    if df_line:
        keywords = [x.lower() for x in df_line.split(';')]
        keywords = list(filter(lambda x: x != "", keywords)) # check non empty strings
        if EC != 0:
            result = 'enzyme'
        elif any(key in class_list for key in keywords):
            for key in keywords:
                if key in class_list:
                    result = key; break
        else:
            result = 'notAnnotated' 
    else:
        result = 'notAnnotated' 
    return result

def get_uniprot_info(protein_list):
    idfrom = "ACC"
    url = 'https://www.uniprot.org/uploadlists/'
    query = " ".join(protein_list)
    #
    ASKCOLUMNS = [  
                    'id',
                    'genes(PREFERRED)',
                    'length',
                    'keywords',
                    'ec'
                    #'families',
                    #'go',
                    #'go-id',
                    #'go(molecular function)'
                    #'go(biological process)',
                    #'go(cellular component)'
                ]
    KEYS = [  
                    'id',
                    'geneid',
                    'length',
                    'keywords',
                    'ec'
                    #'families',
                    #'go',
                    #'go-id',
                    #'go_molfunc'
                    #'go(biological process)',
                    #'go(cellular component)'
                ]
    #
    params = {
                'from': idfrom, 
                'to': 'ACC',
                'format': 'tab',
                'query': query,
                'columns':','.join(ASKCOLUMNS),
                }
    #
    r = requests.post(url, data=params)
    reader = csv.reader(r.text.splitlines(), delimiter = "\t" )
    next(reader)
    #
    records = [
                {key: val for key,val in zip(KEYS, line)}
                for line in reader
                ]
    df_uniprot = pd.DataFrame.from_dict(records).replace('', np.nan)
    return df_uniprot

#####

# 

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)


##########################
########## CODE ##########
OUT_PATH = os.path.join(os.getcwd(), '../Data/pkls')

# List of all proteins
list_proteins = np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()

# list of proteins w structure
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins_struct = list_pdbs + list_afold
logging.info(f'Will retrieve information for {len(list_proteins)}')

class_list = np.loadtxt('../Data/keys.txt', dtype='str', delimiter='\n').tolist()
class_list = [key.lower() for key in class_list]
logging.debug("""keywords from: https://www.uniprot.org/keywords/KW-9992; 
                deleting enzyme keywords because using EC Code for those""")

###############
# Retrieve info from Uniprot
df_uniprot = get_uniprot_info(list_proteins)

df_uniprot.ec = df_uniprot.ec.fillna(0)
df_uniprot.keywords = df_uniprot.keywords.fillna(0)
df_uniprot['Class'] = df_uniprot.apply(lambda x: get_key(x['keywords'], x['ec'], class_list), axis=1)
# clear enzyme
ec_series = df_uniprot.ec[df_uniprot.ec != 0].apply(lambda x: x[:3])
df_uniprot.ec = ec_series

## Source
dic_alfa = dict(zip(list_afold, ['AFold'] * len(list_afold)))
dict_pdb = dict(zip(list_pdbs, ['PDB'] * len(list_pdbs)))
dic_source = {}
dic_source.update(dic_alfa)
dic_source.update(dict_pdb)
df_uniprot['source'] = df_uniprot.id.map(dic_source)

df_uniprot = df_uniprot.drop(columns = ['geneid', 'keywords'])

#
frequencies = df_uniprot['Class'].value_counts()
condition = frequencies < 10
mask_obs = frequencies[condition].index
mask_dict = dict.fromkeys(mask_obs, 'other')
df_uniprot['Class'] = df_uniprot['Class'].replace(mask_dict) 

# counts for cristal
df_uniprot_cristal = df_uniprot[~df_uniprot.source.isna()]
df_uniprot_cristal['Class'].value_counts()

# Save pickles
df_uniprot.to_pickle('../Data/pkls/proteins_annot_full.pkl')
df_uniprot_cristal.to_pickle('../Data/pkls/proteins_annot_crys.pkl')

# df_uniprot_cristal = pd.read_pickle('../Data/pkls/proteins_annot_crys.pkl')
# import matplotlib.pyplot as plt
# plt.clf()
# fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True)
# df_uniprot['Class'].value_counts().plot.bar(ax=axes[0])
# df_uniprot_cristal['Class'].value_counts().plot.bar(ax=axes[1])
# plt.savefig('test_distr_prot.png')
