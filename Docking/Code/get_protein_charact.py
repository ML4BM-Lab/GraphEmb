import os
import numpy as np
import pandas as pd
import glob
import logging
import requests
from tqdm import tqdm
import re 

logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)


# Load objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins = list_pdbs + list_afold

logging.info(f'Retrieving information for {len(list_proteins)}')






data_uniprot = pd.DataFrame(columns=['UniprotID', 'seq_len', 'keyword', 'GOs', 'EC'])


for uniprotid in tqdm(list_proteins):
    r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
    a = r.text.split('\n')
    # SQ INFO
    SQ_info = None
    check_SQ = any('SQ' in string for string in a)
    if check_SQ:
        SQ_info = [line.split() for line in a if 'SQ ' in line][0][2]
    # GOs
    GOs = None
    check_GO = any('GO; ' in string for string in a)
    if check_GO:
        GO_info = [line.strip() for line in a if 'DR   GO;' in line]
        GOs = return_GO(GO_info)
        #GOs
    # ECODE ##--> comprobar si está cogiendo el que le corresponde realmente
    ECode = '-'
    check_E = any('EC=' in string for string in a)
    if check_E: 
        E_Code_line = [line.split() for line in a if 'EC=' in line][0]
        ECode = [element for element in E_Code_line if '=' in element]
        ECode = ECode[0].replace('EC=', '').strip(';')
    # KEYWORDS
    check_KW = any('KW' in string for string in a)
    if check_KW:
        KW_info = [line.replace('KW', '').replace('.', '').split(';') for line in a if 'KW  ' in line]
        keywords = [x.strip().lower() for list_ in KW_info for x in list_]
        keywords = list(filter(lambda x: x != "", keywords)) # check non empty strings
        keywords = [re.sub('.\{(.*?)\}.*', '', key, flags=re.DOTALL) for key in keywords] # remove weird codes
        # special cases 
        if 'ion channel' in keywords: 
            keywords = 'ion channel'
        if 'g-protein coupled receptor' in keywords:
            keywords = 'gpcr'
        if 'motor protein' in keywords:
            keywords = 'motor protein'
    data_uniprot.loc[len(data_uniprot.index)] = [uniprotid, gene_name, SQ_info, keywords, GOs, ECode]


def return_GO(GO_info):
    golist = []
    for line in GO_info:
        gocode = line.split(';')[1].strip()
        gotext = line.split(';')[2].strip().split(':')[1]
        if 'nuclear receptor' in gotext:
            save = 'NR'
            break
        elif 'chaperone binding' in gotext:
            save = 'chaperone binding'
            break
        else:
            golist.append((gocode, gotext))
            #save = golist
            save = None
    return save








df = data_uniprot

df.EC = df.EC.apply(lambda x: str(x)[0])
#df.keyword[df.EC != '-'] = 'Enzyme' # change all that have Ez to Enzyme

# flat_keys = [x for xs in df.keyword[df.keyword != 'Enzyme'].tolist() for x in xs]
# flat_keys =  set(flat_keys) - set(pop_words)

# specify procedency PDB/A fold
df['Source'] = '-'

# creating a Type os Protein column
df['PType'] = np.nan
df['PType'][df.EC != '-']  = 'Enzyme'
df['PType'][df.GOs == 'NR'] = df.GOs[df.GOs == 'NR']
df['PType'][df.GOs == 'chaperone binding'] = df.GOs[df.GOs == 'chaperone binding']
df['PType'][df.keyword == 'gpcr'] = df.keyword[df.keyword == 'gpcr']
df['PType'][df.keyword == 'ion channel'] = df.keyword[df.keyword == 'ion channel']
df['PType'][df.keyword == 'motor protein'] = df.keyword[df.keyword == 'motor protein']


df.PType.shape[0] - df.PType.dropna().shape[0] 

logging.info(f"# Enzymes: {df[df.keyword == 'Enzyme'].shape[0]}")
logging.info(f"# GPCR: {df[df.keyword == 'gpcr'].shape[0]}")
logging.info(f"# IC: {df[df.keyword == 'ion channel'].shape[0]}")
now = df[df.keyword == 'Enzyme'].shape[0] - df[df.keyword == 'gpcr'].shape[0] - df[df.keyword == 'ion channel'].shape[0]
logging.info(f'unclassified prot {df.shape[0] - now }')




