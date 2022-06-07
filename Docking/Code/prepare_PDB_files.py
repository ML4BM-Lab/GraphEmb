
import os
import glob
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from tqdm import tqdm 
import logging


logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)

def clean_pdb(IN_PATH, OUT_PATH, uni, data):
    # use biopython for cleanning
    pdbid = data.loc[uni]['PDB']
    parser = PDBParser()
    path_file = os.path.join(IN_PATH, pdbid+'.pdb')
    structure = parser.get_structure(pdbid, path_file)
    # dip into the structure
    model = structure[0] # for X-Ray is first! 
    chain_id = data.loc[uni]['chain']
    # seg_in =  int(float(data.loc[pdbid]['seg'][0]))
    # seg_out =  int(float(data.loc[pdbid]['seg'][1]))
    chain = model[chain_id]
    # select only the peptidic chain! 
    for residue in list(chain):
        if residue.id[0] != ' ': 
            chain.detach_child(residue.id)
    io = PDBIO()
    io.set_structure(chain)
    io.save(os.path.join(OUT_PATH, f'{uni}.pdb'))



# this new version should take into account
# using uniprot as index
# check if exist the needed pdb

# avoid biopython warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


# Prepare PDB files
IN_PATH = '../Data/PDB_FILES/'
OUT_PATH = "../Data/Clean_from_PDB/"

# available pdbs from pdb
# gunzip * 
# list of pdb files to prepare 
pdb_paths = glob.glob(os.path.join(IN_PATH, '*.pdb'))
pdb_ids = [file.replace(IN_PATH, '').replace('.pdb', '') for file in pdb_paths]

list_download_pdb = np.loadtxt('../Data/list_download_pdbs.txt', dtype=str, delimiter=',').tolist()
len(pdb_ids)

# extract data from chains etcc
res_uni2pdb = pd.read_pickle('../Data/res_uni2pdb.pkl')
res_uni2pdb[['chain', 'seg']] = res_uni2pdb['chain-resid'].str.split('=', n=1, expand=True)
res_uni2pdb['chain'] = res_uni2pdb.chain.apply(lambda x: x[0])
res_uni2pdb['seg'] = res_uni2pdb.seg.str.split('-', n=1, expand=False)
res_uni2pdb = res_uni2pdb.drop(columns = 'chain-resid')

data = res_uni2pdb
data = data.set_index('UniprotID')
for uni in tqdm(data.index):
    if data.loc[uni]['PDB'] in pdb_ids:#print(uni)
        #print('do stuff here')
        clean_pdb(IN_PATH, OUT_PATH, uni, data)
        # clean with the given data
    else:
        logging.debug(f'pdb for uni {uni} not found')
        #continue 



# also need a final list with all retrieved uniprots! 
# with that go to alphafols
logging.debug('0')