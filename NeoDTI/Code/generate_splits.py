import pandas as pd
import os
import sys
sys.path.append('DDR/Code/')
from Model_Sp_Sd_St_split_Improved import generate_splits
import pickle
import numpy as np


#Lets use Yamanishi NR as an example
##Load dataset
#fpath = os.path.join(os.getcwd(),'DB','Data','Yamanashi_et_al_GoldStandard','IC','interactions','ic_admat_dgc_mat_2_line.txt')
#fpath = os.path.join(os.getcwd(),'DB','Data','Yamanashi_et_al_GoldStandard','NR','interactions','nr_admat_dgc_mat_2_line.txt')
#fpath = os.path.join(os.getcwd(),'DB','Data','Davis_et_al','tdc_package_preprocessing','DAVIS_et_al_2line.tsv')
#fpath = os.path.join(os.getcwd(), 'DB', 'Data', 'BIOSNAP', 'ChG-Miner_miner-chem-gene', 'ChG-Miner_miner-chem-gene.tsv')

def mat_to_2line(dti):
    #
    npmat = np.array(dti)
    drugs = [f'D{i}' for i in range(npmat.shape[0])]
    prots = [f'P{i}' for i in range(npmat.shape[1])]
    dp_2line = []
    #
    for i in range(npmat.shape[0]):
        for j in range(npmat.shape[1]):
            if npmat[i,j]:
                dp_2line.append([drugs[i],prots[j]])
    #
    return pd.DataFrame(dp_2line)

##NEODTI (we need to reconstruct the DTI from the drug_protein.txt file )

fpath = os.path.join(os.getcwd(),'NeoDTI','Data','Yamanashi_et_al_GoldStandard','GPCR','mat_drug_protein.txt')
DTIs = mat_to_2line(pd.read_csv(fpath, sep='\t', header = None)) ## MAKE SURE THE HEADER OPTION IS ON/OFF DEPENDING ON THE DATASET!
DTIs.columns = ['Drug', 'Protein']

st_splits = generate_splits(DTIs, mode= 'St', subsampling=True, foldnum=10, 
                            negative_to_positive_ratio = 10, cvopt=True, cvopt_no_names=True,
                            RMSD_dict_opt=False, include_diagonal_RMSD=False)

from collections import Counter

for seed in sp_splits:
    for fold in seed:
        for split in fold:
            a = Counter(list(map(lambda x: x[-1], split)))
            if a[0] == 0 or a[1] == 0:
                print("Imbalanced!")
            
        
    

with open('/mnt/md0/data/jfuente/DTI/Input4Models/NeoDTI/Data/Yamanashi_et_al_GoldStandard/GPCR/gpcr_splits_St.pickle', 'wb') as handle:
    pickle.dump(st_splits, handle, protocol=2)