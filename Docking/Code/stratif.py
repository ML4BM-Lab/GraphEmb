import os
import numpy as np
import pandas as pd
import logging 
import glob

from torch import threshold

data_rmsd = pd.read_csv('../Data/matrix_rmsd_NR.tsv', sep="\t")
list_prot = [colname.replace('_A', '') for colname in data_rmsd.columns.tolist()]
data_rmsd.columns, data_rmsd.index = list_prot, list_prot

data_rmsd

DTIs = pd.read_csv('../../DTINet/Data/Yamanashi_et_al_GoldStandard/NR/final_dtis_NR.tsv', sep = "\t", names=['Drug', 'Protein'])

# using similar to model sp sd...

def get_interactions(DTIs):
    Proteins = DTIs.Protein.unique().tolist()
    Drugs = DTIs.Drug.unique().tolist()
    interactions = []
    for drug in Drugs:
        drug_list = []
        drug_list.append(drug)
        # positive links
        positives =  DTIs.Protein[DTIs.Drug == drug].tolist()
        drug_list.append(positives)
        negatives = list(set(Proteins).difference(set(positives)))
        drug_list.append(negatives)
        # finish list
        interactions.append(drug_list)
    return interactions


interactions = get_interactions(DTIs)

prots_yam = DTIs.Protein.unique().tolist()

#######
prot_pdb = glob.glob('../Data/Clean_from_PDB/*.pdb')
prot_pdb = [prot.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for prot in prot_pdb]
prot_alpha = glob.glob('../Data/Clean_from_AFold/*.pdb')
prot_alpha = [prot.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for prot in prot_alpha]
aval_prots = prot_pdb + prot_alpha
prots_w_str = set(prots_yam).intersection(set(aval_prots))
logging.info(f'Avaliable proteins with structure: {len(prots_w_str)}')

# yamanishi same type of proteins
# let us take a lower threshold zB 2
threshold1 = 3
mask_data = data_rmsd < threshold1

# for only one interacting protein
case = 0
interactions[3][1]
drug = case[0]
positives = case[1]
negatives = case[2]

DTIs_structural = DTIs[DTIs.Protein.isin(prots_w_str)]
interactions = get_interactions(DTIs_structural)


interactions[0]
positive = interactions[0][1]
negatives = interactions[0][2]

result_ = data_rmsd[mask_data[positive]][positive].dropna().drop(positive).sort_values(by=positive).reset_index()
result_.columns = ['uniprotid', 'rmsd']
result_