
import os
import glob
import pandas as pd
import numpy as np


# this new version should take into account
# using uniprot as index
# check if exist the needed pdb



# Prepare PDB files
PATH = '../Data/PDB_FILES/'

# available pdbs from pdb
# list of pdb files to prepare 
pdb_paths = glob.glob(os.path.join(PATH, '*.pdb'))
pdb_ids = [file.replace(PATH, '').replace('.pdb', '') for file in pdb_paths]

len(pdb_ids)

# use biopython for cleanning
res_uni2pdb = pd.read_pickle('../Data/res_uni2pdb.pkl')
res_uni2pdb[['chain', 'seg']] = res_uni2pdb['chain-resid'].str.split('=', n=1, expand=True)
res_uni2pdb['chain'] = res_uni2pdb.chain.apply(lambda x: x[0])
res_uni2pdb['seg'] = res_uni2pdb.seg.str.split('-', n=1, expand=False)
res_uni2pdb = res_uni2pdb.drop(columns = 'chain-resid')

data = res_uni2pdb

data = data.set_index('UniprotID')

for uni in data.index[:10]:
    #print(uni)
    data.loc[uni]
    # check if it was actually downloaded
    if data.loc[uni]['PDB'] in pdb_ids:
        print('do your shite here')
        # function similar to clean pdb
        # but data structure changed ! 
    else:
        print(f'pdb for uni {uni} not found')
        continue


# also need a final list with all retrieved uniprots! 
# with that go to alphafols
