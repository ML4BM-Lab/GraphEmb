import os
import glob
import pandas as pd
import numpy as np

# Prepare PDB files
PATH = '../Data/PDB_FILES/'

###
# 
# https://www.biostars.org/p/251583/ # warnings

# import warnings
# from Bio import BiopythonWarning
# warnings.simplefilter('ignore', BiopythonWarning)

# list of pdb files to prepare 
# glob
pdb_paths = glob.glob(os.path.join(PATH, '*.pdb'))
pdb_ids = [file.replace(PATH, '').replace('.pdb', '') for file in pdb_paths]

len(pdb_ids)
# use biopython for cleanning
res_uni2pdb = pd.read_pickle('../Data/res_uni2pdb.pkl')
res_uni2pdb = res_uni2pdb.set_index('PDB')

list_download = np.loadtxt('../Data/list_download_pdbs.txt', dtype='str', delimiter=',').tolist()
# re download correct dataset! 


res_uni2pdb.loc['1JYK']['chain-resid']
res_uni2pdb[['chain', 'seg']] = res_uni2pdb['chain-resid'].str.split('=', n=1, expand=True)
res_uni2pdb['chain'] = res_uni2pdb.chain.apply(lambda x: x[0])
res_uni2pdb['seg'] = res_uni2pdb.seg.str.split('-', n=1, expand=False)
res_uni2pdb = res_uni2pdb.drop(columns = 'chain-resid')

# floaten segs for selecting segments

# 5LE5
# loop should be over selected uniprots & make them load the selected PDB
# with the selected chain
# sligtly different code

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from tqdm import tqdm 

parser = PDBParser()

OUT_PATH = "../Data/Clean_from_PDB/"

def clean_pdb(IN_PATH, OUT_PATH, pdbid, data):
    # give the path file
    path_file = os.path.join(IN_PATH, pdbid+'.pdb')
    # parse the structure
    structure = parser.get_structure(pdb_ids[0], path_file)
    # dip into the structure
    model = structure[0] # for X-Ray is first! 
    chain_id = data.loc[pdbid]['chain']
    # seg_in =  int(float(data.loc[pdbid]['seg'][0]))
    # seg_out =  int(float(data.loc[pdbid]['seg'][1]))
    chain = model[chain_id]
    # select only the peptidic chain! 
    for residue in list(chain):
        if residue.id[0] != ' ': 
            chain.detach_child(residue.id)
    io = PDBIO()
    io.set_structure(chain)
    io.save(os.path.join(OUT_PATH, f'{pdbid}.pdb'))

for pdbid in tqdm(res_uni2pdb.index):
    #print(pdbid)
    clean_pdb(PATH, OUT_PATH, pdbid, res_uni2pdb )


# if raise waring write it with the name of the pdb? 
# warnings for discontinuous chain ? ? 
