import os
import numpy as np
import pandas as pd
import glob
import logging
import requests
from tqdm import tqdm
from sqlalchemy import create_engine
from sqlalchemy import inspect
import pandas as pd
from sqlalchemy.sql import text
import os
import logging

#####

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)


##########################
########## CODE ##########
# Load objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins = list_pdbs + list_afold

logging.info(f'Retrieving information for {len(list_proteins)}')

class_list = ['ion channel', 'g-protein coupled receptor', 'hormone',
                'cytokine', 'developmental protein', 'chaperone', 'activator'
                'antiviral protein', 'receptor', 'transport', 'transcription']

class_list = ['ion channel', 'g-protein coupled receptor', 'hormone', 'growth factor',
                'motor protein', 
            'transport', 'chaperone', 'transcription']

data_uniprot = pd.DataFrame(columns=['UniprotID', 'seq_len', 'Class', 'EC'])

for uniprotid in tqdm(list_proteins, desc='Looking up in Uniprot'):
    r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
    a = r.text.split('\n')
    SQ_info = None # sequence info
    check_SQ = any('SQ' in string for string in a)
    SQ_info = [line.split() for line in a if 'SQ ' in line if check_SQ][0][2]
    ECode = None # enzyme code info
    check_E = any('EC=' in string for string in a)
    check_BR = any('BRENDA;' in string for string in a)
    if check_E: 
        E_Code_line = [line.split() for line in a if 'EC=' in line][0]
        ECode = [element for element in E_Code_line if '=' in element]
        ECode = ECode[0].replace('EC=', '').strip(';')
    elif check_BR:
        BR_Code_line = [line.split() for line in a if 'BRENDA;' in line]
        ECode = BR_Code_line[0][2].strip(';')
    check_KW = any('KW' in string for string in a) # find class
    if check_KW:
        KW_info = [line.replace('KW', '').replace('.', '').split(';') for line in a if 'KW  ' in line]
        keywords = [x.strip().lower() for list_ in KW_info for x in list_]
        keywords = list(filter(lambda x: x != "", keywords)) # check non empty strings
        if ECode:
            key = 'enzyme'
        elif any(key in class_list for key in keywords):
            for key in keywords:
                if key in class_list:
                    result = key
                    break
        else:
            key = 'notAnnotated' 
            print(uniprotid, keywords, end='\n')
            print('\n')
    data_uniprot.loc[len(data_uniprot.index)] = [uniprotid, SQ_info, key, ECode]


# https://www.uniprot.org/keywords/9999 # opcon biological process

data_uniprot.to_pickle('../Data/data_uniprot_raw.pkl')

# mejor molecular function
# options

# P14649 miosin or motor protein? ['3d-structure', 'motor protein', 'muscle protein', 'myosin', 'reference proteome', 'repeat']
# P0CG47 nucleus ? ['3d-structure', 'adp-ribosylation', 'cytoplasm', 'direct protein sequencing', 'isopeptide bond', 'membrane', 'mitochondrion', 'mitochondrion outer membrane', 'nucleus', 'phosphoprotein', 'reference proteome', 'repeat', 'ubl conjugation']
# add membrane as last struff
# or transmembrane before
# P08648 receptor? ['3d-structure', 'calcium', 'cell adhesion', 'cell junction', 'cleavage on pair of basic residues', 'direct protein sequencing', 'disulfide bond', 'glycoprotein', 'host cell receptor for virus entry', 'host-virus interaction', 'integrin', 'membrane', 'metal-binding', 'phosphoprotein', 'receptor', 'reference proteome', 'repeat', 'signal', 'transmembrane', 'transmembrane helix']
# signal or glicoprotein


## play stuff
# import matplotlib.pyplot as plt

# data = data_uniprot
# https://www.uniprot.org/keywords/KW-9992 igual deberia usar algo d esto
plt.clf()
data_uniprot['Class'].value_counts().plot(kind='pie')
plt.savefig('../Results/test_distr_test2.png')

# data.EC.fillna(value=np.nan, inplace=True)

# # select only first class
# data.EC  = data.EC.dropna().apply(lambda x: str(x)[0]).astype(int)
# #
# data = data.sort_values(by='Class').reset_index()

# mat_or = pd.DataFrame(matrix, columns=ordered, index=ordered)

# data_uniprot.EC.map(enz_dic)

