import os
import glob
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from tqdm import tqdm 
import logging
import warnings
from Bio import BiopythonWarning
import subprocess as sp



def clean_pdb(IN_PATH, OUT_PATH, uni, data):
    pdbid = data.loc[uni]['PDB']
    parser = PDBParser()
    path_file = os.path.join(IN_PATH, pdbid+'.pdb')
    structure = parser.get_structure(pdbid, path_file)
    model = structure[0] # for X-Ray only one model 
    chain_id = data.loc[uni]['chain']
    chain = model[chain_id]
    # select only residues in the peptidic chain
    for residue in list(chain):
        if residue.id[0] != ' ': 
            chain.detach_child(residue.id)
    # safety check 
    if len(list(chain)) > 0:  
        io = PDBIO()
        io.set_structure(chain)
        io.save(os.path.join(OUT_PATH, f'{uni}.pdb'))



##########################################
############### START CODE ###############

# this new version should take into account
# using uniprot as index
# check if exist the needed pdb
# avoid biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)
logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)


# Set Paths
IN_PATH = '../Data/PDB_FILES'
OUT_PATH = "../Data/Clean_from_PDB"

# Unzip pdbs
logging.info('unziping files')
try:
    return_code = sp.check_call(f'gunzip ./{IN_PATH}/*.pdb.gz', shell=True)
    if return_code ==0: 
        logging.info('EXIT CODE 0')
except sp.CalledProcessError as e:
    logging.info(e.output)


# list of pdb files to prepare 
pdb_paths = glob.glob(os.path.join(IN_PATH, '*.pdb'))
pdb_ids = [file.replace(f'{IN_PATH}/', '').replace('.pdb', '') for file in pdb_paths]

list_download_pdb = np.loadtxt('../Data/list_download_pdbs.txt', dtype=str, delimiter=',').tolist()
logging.debug(f'len pdb_ids: {len(pdb_ids)}')

# Extract chain information
res_uni2pdb = pd.read_pickle('../Data/pkls/res_uni2pdb.pkl')
res_uni2pdb[['chain', 'seg']] = res_uni2pdb['chain-resid'].str.split('=', n=1, expand=True)
res_uni2pdb['chain'] = res_uni2pdb.chain.apply(lambda x: x[0])
res_uni2pdb['seg'] = res_uni2pdb.seg.str.split('-', n=1, expand=False)
res_uni2pdb = res_uni2pdb.drop(columns = 'chain-resid')


data = res_uni2pdb
data = data.set_index('UniprotID')
for uni in tqdm(data.index):
    if data.loc[uni]['PDB'] in pdb_ids:
        clean_pdb(IN_PATH, OUT_PATH, uni, data)
    else:
        logging.debug(f'pdb file for uni {uni} not found')



# also need a final list with all retrieved uniprots! 
# with that go to alphafols
logging.debug('0')