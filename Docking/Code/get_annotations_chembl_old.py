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


##
def get_dict_uni2chembl():
    url_mapping = 'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_30/chembl_uniprot_mapping.txt'
    r = requests.get(url_mapping)
    a = r.text.split('\n')
    data = [tuple(cont.split('\t')[:2]) for cont in a[1:-1]] 
    uni2chembl  = dict(data)
    return uni2chembl

#####

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)


# Load objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins = list_pdbs + list_afold

logging.info(f'Retrieving information for {len(list_proteins)}')

logging.debug('Reading chembl_30.db ...')
PATH_DB = '../Data/chembl_30/chembl_30_sqlite'
DB_FILE = 'chembl_30.db'

engine = create_engine(f'sqlite:///{os.path.join(PATH_DB, DB_FILE)}')

# query 
# select all possible levels of classification 
# latter prbibaly  using the first one
qtext = """
SELECT
    target_dictionary.chembl_id          AS target_chembl_id,
    protein_family_classification.l1     AS l1,
    protein_family_classification.l2     AS l2,
    protein_family_classification.l3     AS l3,
    protein_family_classification.l4     AS l4,
    protein_family_classification.l5     AS l5,
    protein_family_classification.l6     AS l6,
    protein_family_classification.l7     AS l7,
    protein_family_classification.l8     AS l8
FROM target_dictionary
    JOIN target_components ON target_dictionary.tid = target_components.tid
    JOIN component_class ON target_components.component_id = component_class.component_id
    JOIN protein_family_classification ON component_class.protein_class_id = protein_family_classification.protein_class_id
"""

with engine.begin() as conn:
    res = conn.execute(text(qtext))
    df = pd.DataFrame(res.fetchall())

df.columns = res.keys()
df = df.where((pd.notnull(df)), None)

logging.debug(f"check nAChr {df[df['target_chembl_id'] == 'CHEMBL4808']}")

# Annotating with first level
logging.info('We are annotating with the first level for family (most general)')

data = df[['target_chembl_id', 'l1']]
uni2chembl = get_dict_uni2chembl()
chembl2uni = dict(zip(uni2chembl.values(), uni2chembl.keys()))

logging.debug('playing with data...')

data.target_chembl_id.map(chembl2uni)
df = pd.DataFrame(list_proteins,  columns=['UniprotID'])
chembl2family = dict(zip(data.target_chembl_id.tolist(), data.l1.tolist()))
df['CHEMBL'] = df.UniprotID.map(uni2chembl)
df['Family'] = df.CHEMBL.map(chembl2family)
df = df.drop(columns='CHEMBL')
df = df.fillna(value = 'NotAnnotated')

## Source
dic_alfa = dict(zip(list_afold, ['AFold'] * len(list_afold)))
dict_pdb = dict(zip(list_pdbs, ['PDB'] * len(list_pdbs)))
dic_source = {}
dic_source.update(dic_alfa)
dic_source.update(dict_pdb)
df['Source'] = df.UniprotID.map(dic_source)


# try with the highest category in  remaing in uniprot
# at the same time retrieve sequene

list_remaining = df[df.Family == 'NotAnnotated'].UniprotID.tolist()


data_uniprot = pd.DataFrame(columns=['UniprotID', 'seq_len', 'keyword', 'EC'])
pop_words = ['3d-structure', 'calcium', 'zinc',
                'mental retardation', 'neurogenesis']
words = []
for uniprotid in tqdm(list_proteins):
    #uniprotid
    r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
    a = r.text.split('\n')
    # SQ INFO
    SQ_info = None
    check_SQ = any('SQ' in string for string in a)
    SQ_info = [line.split() for line in a if 'SQ ' in line if check_SQ][0][2]
    ECode = None
    check_E = any('EC=' in string for string in a)
    check_BR = any('BRENDA;' in string for string in a)
    if check_E: 
        E_Code_line = [line.split() for line in a if 'EC=' in line][0]
        ECode = [element for element in E_Code_line if '=' in element]
        ECode = ECode[0].replace('EC=', '').strip(';')
    elif check_BR:
        BR_Code_line = [line.split() for line in a if 'BRENDA;' in line]
        ECode = BR_Code_line[0][2].strip(';')
        #Code
    # KEYWORDS
    check_KW = any('KW' in string for string in a)
    if check_KW:
        KW_info = [line.replace('KW', '').replace('.', '').split(';') for line in a if 'KW  ' in line]
        keywords = [x.strip().lower() for list_ in KW_info for x in list_]
        keywords = list(filter(lambda x: x != "", keywords)) # check non empty strings
        for key in keywords:
            if key in pop_words: keywords.remove(key)
            # special cases 
        if ECode:
            keyword = 'enzyme'
        elif 'ion channel' in keywords: 
            keyword = 'ion channel'
        elif 'g-protein coupled receptor' in keywords:
            keyword = 'gpcr'
        elif 'transport' in keywords:
            keyword = 'transport'
        elif 'chaperone' in keywords:
            keyword = 'chaperone'
        elif 'hormone' in keywords:
            keyword = 'hormone'
        elif 'transcription' in keyword:
            keyword = 'transcription'
        else:
            keyword = 'notAnnotated' #keywords
            words.append(keywords)
        #print('\n')
        #time.sleep(3)
    data_uniprot.loc[len(data_uniprot.index)] = [uniprotid, SQ_info, keyword, ECode]




data_uniprot[data_uniprot.keyword != 'Enzyme'].iloc[0]

data_uniprot[data_uniprot.keyword != 'Enzyme'].iloc[2].keyword

import matplotlib.pyplot as plt
plt.clf()
data_uniprot['keyword'].value_counts().plot(kind='pie')
plt.savefig('test_distr.png')


data_uniprot.keyword[data_uniprot.keyword != 'Enzyme']

for i, j in zip(data_uniprot.UniprotID, data_uniprot.keyword):
    i,j



test_ = data_uniprot