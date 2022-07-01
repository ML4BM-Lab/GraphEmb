import os
import numpy as np
import pandas as pd
import glob
import logging
import requests
from tqdm import tqdm
import re 
from sqlalchemy import create_engine
from sqlalchemy import inspect
import pandas as pd
from sqlalchemy.sql import text
import os
import logging
import csv
import time


##
def get_dict_chembl2uni():
    logging.info('Retrieving dictionary chembl2uni from rest api')
    # original: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_30/chembl_uniprot_mapping.txt
    PATHDICT = '../Data/chembl_uniprot_mapping.txt'
    with open(PATHDICT) as f:
        a = f.readlines()
    a = [line.strip('\n') for line in a]
    data = [tuple(cont.split('\t')[:2]) for cont in a[1:]] 
    chembl2uni = {chembl:uni for uni,chembl in data}
    return chembl2uni


def retrieve_uniprot_data(protein_list):
    idfrom = "ACC"
    url = 'https://www.uniprot.org/uploadlists/'
    query = " ".join(protein_list)
    #
    ASKCOLUMNS = [  
                    'id',
                    'entry name',
                    'genes(PREFERRED)',
                    'length',
                    'ec',
                    'families'
                ]
    KEYS = [  
                    'id',
                    'entry name',
                    'genes_preferred',
                    'length',
                    'ec',
                    'families'
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
    #
    df_uniprot = pd.DataFrame.from_dict(records).replace('', np.nan)
    return df_uniprot



#########

logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)


####################################
####################################


# Load objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins_str = list_pdbs + list_afold

list_proteins = np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()

logging.info(f'Retrieving information for {len(list_proteins)}')

#######################################
##### CHEMBL

logging.debug('Reading chembl_30.db ...')
# Paths
PATH_DB = '../Data/chembl_30/chembl_30_sqlite'
DB_FILE = 'chembl_30.db'

# load
engine = create_engine(f'sqlite:///{os.path.join(PATH_DB, DB_FILE)}')

# query: select all possible levels of classification 
# qtext = """
# SELECT
#     target_dictionary.chembl_id          AS target_chembl_id,
#     protein_family_classification.l1     AS l1,
#     protein_family_classification.l2     AS l2,
#     protein_family_classification.l3     AS l3,
#     protein_family_classification.l4     AS l4,
#     protein_family_classification.l5     AS l5,
#     protein_family_classification.l6     AS l6,
#     protein_family_classification.l7     AS l7,
#     protein_family_classification.l8     AS l8
# FROM target_dictionary
#     JOIN target_components ON target_dictionary.tid = target_components.tid
#     JOIN component_class ON target_components.component_id = component_class.component_id
#     JOIN protein_family_classification ON component_class.protein_class_id = protein_family_classification.protein_class_id
# """

qtext = """
SELECT
    target_dictionary.chembl_id          AS target_chembl_id,
    protein_family_classification.l1     AS l1,
    protein_family_classification.l2     AS l2,
    protein_family_classification.l3     AS l3,
    protein_family_classification.l4     AS l4
FROM target_dictionary
    JOIN target_components ON target_dictionary.tid = target_components.tid
    JOIN component_class ON target_components.component_id = component_class.component_id
    JOIN protein_family_classification ON component_class.protein_class_id = protein_family_classification.protein_class_id
"""

with engine.begin() as conn:
    res = conn.execute(text(qtext))
    df_chembl = pd.DataFrame(res.fetchall())

df_chembl.columns = res.keys()
df_chembl = df_chembl.where((pd.notnull(df_chembl)), None)


logging.debug(f"check nAChr {df_chembl[df_chembl['target_chembl_id'] == 'CHEMBL4808']}")

chembl2l1 = dict(zip(df_chembl.target_chembl_id, df_chembl.l1))
chembl2l2 = dict(zip(df_chembl.target_chembl_id, df_chembl.l2))

# Annotating with first level
logging.info('We are annotating with the first level for family (most general)')

# get dict
chembl2uni = get_dict_chembl2uni()

## Source
dic_alfa = dict(zip(list_afold, ['AFold'] * len(list_afold)))
dict_pdb = dict(zip(list_pdbs, ['PDB'] * len(list_pdbs)))
dic_source = {}
dic_source.update(dic_alfa)
dic_source.update(dict_pdb)


##### batch access to uniprot
df_annot_uni = retrieve_uniprot_data(list_proteins)
df_annot_uni

#### create a new db
df = df_annot_uni[['id', 'length', 'ec']].copy()

df['source'] = df.id.map(dic_source)

uni2chembl = dict(zip(chembl2uni.values(), chembl2uni.keys()))
df['chembl'] = df.id.map(uni2chembl)
df['l1'], df['l2'] = df.chembl.map(chembl2l1), df.chembl.map(chembl2l2)

# creating dicctionaries
ec_series = df_annot_uni.ec[~df_annot_uni.ec.isna()].apply(lambda x: x[:3])
df_annot_uni.ec = ec_series


######
# fast vis
import matplotlib.pyplot as plt

df.l1[~df.source.isna()].value_counts().plot(kind='pie')
plt.savefig('test_data.pdf')
plt.clf()

df.l1 = df.l1.fillna('Unclassified protein')

df[~df.source.isna()].l1.isna().sum()

#