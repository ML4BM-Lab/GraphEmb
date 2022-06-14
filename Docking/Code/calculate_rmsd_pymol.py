import glob
from pymol import cmd
from align_all_to_all import align_all_to_all
import time
from tqdm import tqdm
import importlib


# do for all available drugs
# launch with nohup

# modified stugf. 

list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_objects = list_pdbs + list_afold

# load pymol objects
for object in tqdm(list_objects):
    cmd.load(object)

# calculate rmsd
print('Calculating rmsd (from superimposition)...')
a = time.time()
align_all_to_all(object_list=None,selection='name ca',cutoff=2,cycles=5,debug=0,full_matrix=1,method='super')
b = time.time()

t = b-a
print(f"{t/60} minutes")

'''
### 
## 2. using super 
# (Superimpose is recommended by Pymol: 
# super is more robust  than align 
# for proteins with low sequence similarity. ) 
# nuclear receptors as example
threshold = 5 # similarity

import matplotlib.pyplot as plt
import seaborn as sns
m = pd.read_csv('test_E.txt', header=0, sep ="\t+", index_col=0,  engine='python')

# read as pd.read_csv('test.txt', header=0, sep ="\t+", index_col=0)

plt.clf()
sns.heatmap(m<10, cmap='RdYlGn')
#m[m[m<5]>0.5]
plt.savefig('test_result.png')

'''
