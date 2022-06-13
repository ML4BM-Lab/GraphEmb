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

uniprotid = 'O75899'

r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
a = r.text.split('\n')

check_pdb = any('PDB;' in string for string in a)
if check_pdb:
    pdb_info = [line for line in a if ' PDB; ' in line]
    for s in pdb_info:
        info = s.split(';')
        info = [nfo.strip() for nfo in info]
        df_uni2pdb.loc[len(df_uni2pdb.index)] = [uniprotid, info[1], info[2], info[3], info[4]]


check_pdb = any('KW' in string for string in a)
KW_info = [line for line in a if 'KW ' in line]


import requests
from xml.etree import ElementTree

response = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.xml')

tree = ElementTree.fromstring(response.content)

response = requests.get('http://blabla.com')
dict_data = xmltodict.parse(response.content)
# xmltodict

for i in dict_data['uniprot']['entry']['keyword']:
    i['#text']