import os
import pandas as pd
import logging
from tqdm import tqdm
import helper_functions as hf
import requests
from tqdm import tqdm
#import glob 
#from collections import Counter
import numpy as np
#from Bio.PDB import PDBP

uniprotid = 'Q9Y698'

#def get_info_uniprot(list_targets):
uniprotid = 'Q9F0I5'

df_uni2pdb = pd.DataFrame(columns=['UniprotID', 'PDB', 'Method', 'Resolution', 'chain-resid', 'keyword', 'EC'])

for uniprotid in tqdm(list_targets, desc=f'Retrieving info for proteins'):
    r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
    a = r.text.split('\n')
    # Check PDB info
    check_pdb = any('PDB;' in string for string in a)
    if check_pdb:
        # loop inside check here fails
        pdb_info = [line for line in a if ' PDB; ' in line]
        for s in pdb_info:
            info = s.split(';')
            info = [i.strip() for i in info]
    # check keywords
    check_KW = any('KW' in string for string in a)
    if check_KW:
        KW_info = [line.replace('KW', '').replace('.', '').split(';') for line in a if 'KW ' in line]
        keywords = [x.strip().lower() for list_ in KW_info for x in list_]
    # check Enzyme Code
    check_E = any('EC=' in string for string in a)
    if check_E:
        E_Code_line = [line.split() for line in a if 'EC=' in line][0]
        E_Code = [element for element in E_Code_line if 'EC=' in element][0].split('=')[-1]
        if ';' in E_Code: E_Code = E_Code.strip(';')
    else:
        E_Code = '-'
        #E_Code
    df_uni2pdb.loc[len(df_uni2pdb.index)] = [uniprotid, info[1], info[2], info[3], info[4], keywords, E_Code]
#    return df_uni2pdb


import xmltodict
import requests

list_targets = all_prot


df_uni2pdb = pd.DataFrame(columns=['uniprotid', 'gene_id', 'pdbid', 'method', 'resolution', 'chains', 'keywords', 'EC_Code'])

for uniprotid in tqdm(list_targets[:10]):
    uniprotid
    url = f'https://www.uniprot.org/uniprot/{uniprotid}.xml'
    try:   
        r = requests.get(url)
        data = xmltodict.parse(r.content)
    except:
        print(f'error trying to parse  {uniprotid}')    
    data_entry = data['uniprot']['entry']
    if 'PDB' in [i.get('@type',0) for i in data_entry['dbReference']]:
        EC_Code = '-'
        if data_entry['protein'].get('recommendedName') and data_entry['protein'].get('recommendedName', {}).get('ecNumber'):
            ec_entry = data_entry['protein']['recommendedName']['ecNumber']
            if type(ec_entry) == str:
                EC_Code = data_entry['protein']['recommendedName']['ecNumber']
            elif type(ec_entry) == list:
                if ec_entry[0] != dict:
                    EC_Code = ec_entry[0]
                else:
                    EC_Code = [element['#text'] for element in data_entry['protein']['recommendedName']['ecNumber']][0]
            else:
                EC_Code = data_entry['protein']['recommendedName']['ecNumber'].get('#text')
        if data_entry.get('gene') and data_entry.get('gene', {}).get('name'):
            gene_name = data_entry['gene']['name']
            if type(gene_name) is list:
                gene_id = [text.get('#text') for text in data_entry['gene']['name']]
            else:
                gene_id = gene_name.get('#text')
        #keywords
        keywords = None
        if type(data_entry['keyword']) == list:
            keywords = [key['#text'] for key in data_entry['keyword']]
        # pdb information
        for reference in data_entry['dbReference']:
            if reference['@type'] == 'PDB':
                pdbid = reference['@id']
                props = reference['property']
                for prop in props:
                    if prop['@type'] == 'method': 
                        method = prop['@value']
                    elif prop['@type'] == 'resolution':
                        resolution = prop['@value']
                    elif prop['@type'] == 'chains':
                        chains = prop['@value']
                df_uni2pdb.loc[len(df_uni2pdb.index)] = [uniprotid, gene_id, pdbid, method, resolution, chains, keywords, EC_Code]


