import os
import pandas as pd
import logging
from tqdm import tqdm
import glob 
import numpy as np
from Bio.PDB import PDBParser
#from Bio.PDB.MMCIFParser import MMCIFParser
import gzip
import shutil
import subprocess as sp



def get_meanLDDT_for_pdb(file):
    parser = PDBParser()
    structure = parser.get_structure('afold_model', gzip.open(file, "rt"))
    avg_bfactor_ = []
    # SMCRA (Structure/Model/Chain/Residue/Atom) format
    model = structure[0] # # X-Ray generally only have 1 model, while more in NMR
    for chain in model:
        for residue in chain:
            for atom in residue:
                avg_bfactor_.append(atom.get_bfactor()) 
    avg_bfactor = np.mean(avg_bfactor_)
    return avg_bfactor

# def get_meanLDDT_for_cif(file):
#     parser = MMCIFParser()
#     structure = parser.get_structure('afold_model', gzip.open(file, "rt"))
#     avg_bfactor_ = []
#     # SMCRA (Structure/Model/Chain/Residue/Atom) format
#     #for model in structure:
#     model = structure[0] # # X-Ray generally only have 1 model, while more in NMR
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 avg_bfactor_.append(atom.get_bfactor()) # 58.60314146341464
#     avg_bfactor = np.mean(avg_bfactor_)
#     return avg_bfactor



##########################################
############### START CODE ###############
logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)

# alpha fold => no option to download or check from the bulk downloaded data
logging.info('Retrieving afold data...')
listdir_aphold = glob.glob('../Data/AFOLD_FILES/*.pdb.gz')
uniprot_aphold_list = list(set([line.split('-')[1] for line in listdir_aphold]))
logging.info(f'Number of available human proteins in alphafold: {len(uniprot_aphold_list)}')

# Instead of calculating the average for all proteins, 
# we compare for those that we need and do not appear in PDB
all_prot = np.loadtxt('../Data/all_prot.txt', dtype=str).tolist()
pdb_prot = glob.glob(os.path.join('../Data/Clean_from_PDB/', '*.pdb')) # Clean files
pdb_prot = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in pdb_prot]
missing_prot = set(all_prot) - set(pdb_prot) # check

# Intersect those available in alphafold
afold_aval = list(set(uniprot_aphold_list).intersection(missing_prot))
logging.info(f'From {len(missing_prot)} we can rescue {len(afold_aval)} from alpha fold')
# read files
afold_aval_files = [ f'../Data/AFOLD_FILES/AF-{file}-F1-model_v2.pdb.gz' for file in afold_aval]
afold_pdb_data = pd.DataFrame(columns=['UniprotID', 'avgpLDDT'])
for file in tqdm(afold_aval_files, desc='Retrieving mean LDDT'):
    uniid = file.split('-')[1]    
    avg_bfactor = get_meanLDDT_for_pdb(file) # with cif takes hours 
    afold_pdb_data.loc[len(afold_pdb_data.index)] = [uniid, avg_bfactor]

# confident > 70
folder_checkpoint = '../Data/pkls'
afold_pdb_data.avgpLDDT = afold_pdb_data.avgpLDDT.round(2)
afold_pdb_data_70 = afold_pdb_data[afold_pdb_data.avgpLDDT>70]
afold_pdb_data_70.to_pickle(os.path.join(folder_checkpoint, 'afold_data_70.pkl'))

# get the selected uniprots and copy files to our work folder
selected_afold_id = afold_pdb_data_70.UniprotID.unique().tolist()
logging.info(f'Will take {len(selected_afold_id)} proteins from alpha fold')
selected_afold_file = [f'../Data/AFOLD_FILES/AF-{sfile}-F1-model_v2.pdb.gz' for sfile in selected_afold_id]

for id, file in zip(selected_afold_id, selected_afold_file):
    original = file
    target = os.path.join('../Data/Clean_from_AFold', f'{id}.pdb.gz')
    shutil.copyfile(original, target)


# Unzip pdbs
logging.info('unziping files')
try:
    return_code = sp.check_call(f'gunzip ./../Data/Clean_from_AFold/*.pdb.gz', shell=True)
    if return_code ==0: 
        logging.info('EXIT CODE 0')
except sp.CalledProcessError as e:
    logging.info(e.output)
