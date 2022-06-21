import time
from tqdm import tqdm
import logging
import glob
from pymol import cmd
from align_all_to_all import align_all_to_all


##########################################
############### START CODE ###############
logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)

# Load objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_objects = list_pdbs + list_afold

# Load Pymol objects
for object in tqdm(list_objects, desc='Loading pdbs'):
    cmd.load(object)

# Calculate RMSD
# Superimpose is recommended by Pymol: 
# super is more robust  than align 
# for proteins with low sequence similarity. ) 
print('Calculating RMSD (method = super)...')
a = time.time()
align_all_to_all(object_list=None,selection='name ca',cutoff=2,cycles=5,debug=0,full_matrix=1,method='super')
b = time.time()

t = b-a
print(f"Time elapsed: {t/60:.4f} minutes")
