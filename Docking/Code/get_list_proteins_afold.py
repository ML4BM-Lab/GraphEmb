import os
import pandas as pd
import logging
from tqdm import tqdm
import helper_functions as hf
import requests
from tqdm import tqdm
import glob 
from collections import Counter
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import gzip
import shutil


def get_meanLDDT_for_pdb(file):
    parser = PDBParser()
    structure = parser.get_structure('afold_model', gzip.open(file, "rt"))
    avg_bfactor_ = []
    # SMCRA (Structure/Model/Chain/Residue/Atom) format
    #for model in structure:
    model = structure[0] # # X-Ray generally only have 1 model, while more in NMR
    for chain in model:
        for residue in chain:
            for atom in residue:
                avg_bfactor_.append(atom.get_bfactor()) # 58.60314146341464
    avg_bfactor = np.mean(avg_bfactor_)
    return avg_bfactor

def get_meanLDDT_for_cif(file):
    parser = MMCIFParser()
    structure = parser.get_structure('afold_model', gzip.open(file, "rt"))
    avg_bfactor_ = []
    # SMCRA (Structure/Model/Chain/Residue/Atom) format
    #for model in structure:
    model = structure[0] # # X-Ray generally only have 1 model, while more in NMR
    for chain in model:
        for residue in chain:
            for atom in residue:
                avg_bfactor_.append(atom.get_bfactor()) # 58.60314146341464
    avg_bfactor = np.mean(avg_bfactor_)
    return avg_bfactor

##
logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)
##

######### ALPHA FOLD
# alpha fold => no option to download or check from the bulk downloaded data
print('getting afold data')
listdir_aphold = glob.glob('../ALPHA_FILES/*.pdb.gz')
uniprot_aphold_list = list(set([line.split('-')[1] for line in listdir_aphold]))
logging.info(f'Number of available uniprots with alphafold {len(uniprot_aphold_list)}')

##
# instead of calculating the average of all proteins, we compare
# those that we need 
all_prot = np.loadtxt('../Data/all_prot.txt', dtype=str).tolist()
pdb_prot = glob.glob(os.path.join('../Data/Clean_from_PDB/', '*.pdb'))
pdb_prot = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in pdb_prot]
missing_prot = set(all_prot) - set(pdb_prot)

# intersect those available in alphafold
afold_aval = list(set(uniprot_aphold_list).intersection(missing_prot))
logging.info(f'From {len(missing_prot)} we can rescue {len(afold_aval)} from alpha fold')

##
# biopython
afold_aval_files = ['../Data/AFOLD_FILES/' + 'AF-' + file + '-F1-model_v2.pdb.gz' for file in afold_aval]

afold_pdb_data = pd.DataFrame(columns=['UniprotID', 'avgpLDDT'])
for file in tqdm(afold_aval_files):
    uniid = file.split('-')[1]    
    avg_bfactor = get_meanLDDT_for_pdb(file) # with cif takes hours and we are just calculating average score
    afold_pdb_data.loc[len(afold_pdb_data.index)] = [uniid, avg_bfactor]

#afold_pdb_data.to_pickle('../Data/afold_all_data.pkl')

# confident > 70
afold_pdb_data.avgpLDDT = afold_pdb_data.avgpLDDT.round(2)
afold_pdb_data_70 = afold_pdb_data[afold_pdb_data.avgpLDDT>70]
afold_pdb_data_70.to_pickle('../Data/afold_data_70.pkl')
selected_afold_id = afold_pdb_data_70.UniprotID.unique().tolist()

selected_afold_file = ['../Data/AFOLD_FILES/' + 'AF-' + file + '-F1-model_v2.pdb.gz' for file in selected_afold_id]

for id, file in zip(selected_afold_id, selected_afold_file):
    #print(id, file)
    original = file
    target = os.path.join('../Data/Clean_from_AFold', f'{id}.pdb.gz')
    #print(original, target)
    shutil.copyfile(original, target)

# and unzip them
# gunzip * 
