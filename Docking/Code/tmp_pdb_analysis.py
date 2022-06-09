import requests
from tqdm import tqdm
import glob

pdb_files = glob.glob('../Data/PDB_FILES/*.pdb')

pdb_files = [file.replace('.pdb', '').replace('../Data/PDB_FILES/','').lower()  for file in pdb_files]


# this is with pdbe
pdb_incomplex = []
for pdbid in tqdm(pdb_files, desc='looking for complex'):
    r = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdbid}')
    if r.status_code != 200: 
        print('nope')
        continue
    title = r.json()[pdbid][0]['title'].lower()
    if 'complex' in title:
        pdb_incomplex.append(pdbid.upper())


##
pdb_files = [file.replace('.pdb', '').replace('../Data/PDB_FILES/','')  for file in pdb_files]

pdb_incomplex_PDB = []

for pdbid in tqdm(pdb_files, desc='looking for complex'):
    r = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{pdbid}')
    if r.status_code != 200: 
        print(f'nope {pdbid}')
        continue
    title = r.json()['struct']['title'].lower()
    if 'complex' in title: 
        pdb_incomplex_PDB.append(pdbid)

