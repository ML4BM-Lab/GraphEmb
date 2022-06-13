import glob
from pymol import cmd
from align_all_to_all import align_all_to_all
import time
from tqdm import tqdm
import importlib
import numpy as np

#filename1 = '../Data/Clean_from_PDB/P09960.pdb'
#filename2 = '../Data/Clean_from_PDB/P52687.pdb'

# glob.glob('Data_clean')
#nr_in_pdb = ['P37231', 'O75469', 'P35398', 'P06401', 'Q96RI1', 'O95718','P10826']
#nr_in_afold = ['Q14541', 'Q14994']

# list_objects = 
print('testing E') 
e_in_afold = np.loadtxt('../Data/test_E_docking/test_E_docking_fromafold.txt', dtype=str).tolist()
e_in_pdb = np.loadtxt('../Data/test_E_docking/test_E_docking_frompdb.txt', dtype=str).tolist()
list_pdb = ['../Data/Clean_from_PDB/' + pdbid + '.pdb' for pdbid in e_in_pdb]
list_afold = ['../Data/Clean_from_AFold/' + pdbid + '.pdb' for pdbid in e_in_afold]
list_objects = list_pdb + list_afold

# test time 
#list_objects = glob.glob('../Data/Clean_from_PDB/*.pdb')

for object in tqdm(list_objects):
    cmd.load(object)

print('Calculating rmsd (via superimposition)...')
a = time.time()
align_all_to_all(object_list=None,selection='name ca',cutoff=2,cycles=5,debug=0,full_matrix=1,method='super')
b = time.time()

t = b-a
print(f"{t/60} minutes")