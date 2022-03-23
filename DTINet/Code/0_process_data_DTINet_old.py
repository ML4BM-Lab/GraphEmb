import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import multiprocessing as mp
import time
import json

# Transform & do this in: GET_ALL_MATRIX.py
# import all coordinates

wdir = '../Data/DrugBank'

files = ['coordinates_drug_disease.tsv', 'coordinates_drug_se.tsv', 'coordinates_PPI.tsv',
            'coordinates_drug_drug.tsv', 'DrugBank_DTIs.tsv',  'coordinates_protein_disease.tsv']

drug_dis = pd.read_csv(os.path.join(wdir, files[0]), sep='\t')
drug_se = pd.read_csv(os.path.join(wdir, files[1]), sep='\t')
PPI = pd.read_csv(os.path.join(wdir, files[2]), sep='\t', names=['P1', 'P2'])
drug_drug = pd.read_csv(os.path.join(wdir, files[3]), sep='\t', names=['D1', 'D2'])
DTI = pd.read_csv(os.path.join(wdir, files[4]), sep='\t', names=['Drug', 'Protein'])
prot_dis = pd.read_csv(os.path.join(wdir, files[5]), sep='\t', usecols=['UniprotID', 'DiseaseID'])


###############################################################################################
################################## GET FILTERED NODES #########################################

# SELECT DRUG NODES
#all_drugs = np.loadtxt(os.path.join(wdir, 'all_drugs.txt'), dtype=str).tolist()
drugs_linked_to_drugs = set(drug_drug.D1.tolist() + drug_drug.D2.tolist()) # 4418
drugs_linked_to_disease = set(drug_dis.DrugBankID.tolist()) # 2754
drugs_linked_to_sideeffect = set(drug_se.DrugBank_ID.tolist()) # 675
drugs_linked_to_proteins = set(DTI.Drug.tolist()) # 7627
not_isolated_drugs = list(drugs_linked_to_drugs.union(drugs_linked_to_disease, drugs_linked_to_sideeffect, drugs_linked_to_proteins )) ### 9720
not_isolated_drugs.sort()
#### need to have the drug_w_smiles.txt first
list_drug_w_smiles = np.loadtxt(os.path.join(wdir, 'drug_w_smiles.txt'), dtype=str).tolist() ##############################------>>> ** aqui ya estaria! 
len(set(not_isolated_drugs).intersection(set(list_drug_w_smiles))) # 8638


### SELECT PROTEIN NODES
#all_proteins = np.loadtxt(os.path.join(wdir, 'all_protein.txt'), dtype=str).tolist() # 9183
proteins_linked_to_proteins = set(PPI.P1.tolist() + PPI.P2.tolist()) # 9183 
proteins_linked_to_disease = set(prot_dis.UniprotID.tolist()) # 8985
proteins_linked_to_drugs = set(DTI.Protein.tolist()) # 4884
not_isolated_proteins = list(proteins_linked_to_proteins.union(proteins_linked_to_disease, proteins_linked_to_drugs)) # 11941
not_isolated_proteins # not any specific order 


##################### ****************************************************************************************************
## quitar las que no tengan smiles y las que no tengan sim en seq (SMILES DONE; MISSING SEQUENCE!)
## ==> intersec entre not isolated & have_smiles
## not_isolated_drugs.intersect(have_smiles)
## => not_isolated_drugs = not_isolated_drugs_w_smiles
##################### *****************************************************************************************************

# reference as in original data. Save both files as
np.savetxt(os.path.join(wdir ,'drug.txt'), not_isolated_drugs , newline='\n', fmt='%s')
np.savetxt(os.path.join(wdir ,'protein.txt'), not_isolated_proteins , newline='\n', fmt='%s')

###############################################################################################
################################### GET PROCESSED MATRIX ######################################

# for matrix, use
#  not_isolated_drugs & not_isolated_proteins 

######## Drug matrix
# DRUG DISEASE MATRIX
# drug_dis
matrix_drug_dis_ = pd.get_dummies(drug_dis.set_index('DrugBankID')['DiseaseID']).max(level=0)#.reset_index()# no esta ordenado esto
matrix_drug_dis_.columns

len(set(matrix_drug_dis_.columns))
matrix_drug_dis_.shape # --> (2754, 7072)

matrix_drug_dis = pd.DataFrame(matrix_drug_dis_, columns= matrix_drug_dis_.columns, index= not_isolated_drugs)
matrix_drug_dis = matrix_drug_dis.fillna(int(0))
matrix_drug_dis = matrix_drug_dis.astype(int)
matrix_drug_dis.head(2)
matrix_drug_dis.to_csv(os.path.join(wdir, 'mat_drug_dis.txt'), index=False, header=False, sep=" ") 

# DRUG SIDE EFFECT MATRIX
matrix_drug_se_ = pd.get_dummies(drug_se.set_index('DrugBank_ID')['se']).max(level=0)
matrix_drug_se_.columns

len(set(matrix_drug_se_.columns))
matrix_drug_se_.shape # --> (675, 4669)

matrix_drug_se = pd.DataFrame(matrix_drug_se_, columns= matrix_drug_se_.columns, index= not_isolated_drugs)
matrix_drug_se = matrix_drug_se.fillna(int(0))
matrix_drug_se = matrix_drug_se.astype(int)
matrix_drug_se.head(2)
matrix_drug_se.to_csv(os.path.join(wdir, 'mat_drug_se.txt'), index=False, header=False, sep=" ") 


# DRUG DRUG MATRIX
matrix_drug_drug_ = pd.get_dummies(drug_drug.set_index('D1')['D2']).max(level=0)
matrix_drug_drug_.columns

len(set(matrix_drug_drug_.columns))
matrix_drug_drug_.shape # --> (4417, 4418)

matrix_drug_drug = pd.DataFrame(matrix_drug_drug_, columns= not_isolated_drugs, index= not_isolated_drugs)
matrix_drug_drug = matrix_drug_drug.fillna(int(0))
matrix_drug_drug = matrix_drug_drug.astype(int)
matrix_drug_drug.head(2)
matrix_drug_drug.to_csv(os.path.join(wdir, 'mat_drug_drug.txt'), index=False, header=False, sep=" ") 



matrix_drug_drug['DB00006']['DB06605']
matrix_drug_drug['DB06605']['DB00006']

####### Protein matrix#################################################################################################### ******** ####################################
###### esta matrix ten o problema de que só estan recollidas as interaccións dun lado, e se a construímos así non vai funcionar 
# é coma se estivese dirixida. X interacciona con Y, pero Y non interacciona con X (ollo)
# entaõ facer P1 - P2 & P2 - P1, e facer fusionar !!! 
# PROTEIN PROTEIN MATRIX ##########################################################################################
matrix_protein_protein_ = pd.get_dummies(PPI.set_index('P1')['P2']).max(level=0)
matrix_protein_protein_.columns

len(set(matrix_protein_protein_.columns))
matrix_protein_protein_.shape # --> (7618, 7166)
matrix_protein_protein = pd.DataFrame(matrix_protein_protein_, columns= not_isolated_proteins, index= not_isolated_proteins)
matrix_protein_protein = matrix_protein_protein.fillna(int(0))
matrix_protein_protein = matrix_protein_protein.astype(int)
matrix_protein_protein.head(2)
matrix_protein_protein.to_csv(os.path.join(wdir, 'mat_protein_protein.txt'), index=False, header=False, sep=" ") 

#  Q13683  P02708
matrix_protein_protein['P02708']['Q13683']

mat_protein_protein_check_columns = matrix_protein_protein.replace(0,np.nan)
ocol = mat_protein_protein_check_columns.shape[1]
ncol = mat_protein_protein_check_columns.dropna(axis=1, how='all').shape[1]

mat_protein_protein_check_rows = matrix_protein_protein.replace(0,np.nan)
orows = mat_protein_protein_check_rows.shape[0]
nrows = mat_protein_protein_check_rows.dropna(axis=0, how='all').shape[0]
print('Drug Protein')
print(f'Lost {ocol-ncol} empty columns (protein)')
print(f'Lost {orows-nrows} empty rows (protein)')


# PROTEIN DISEASE MATRIX

prot_dis
matrix_prot_dis_ = pd.get_dummies(prot_dis.set_index('UniprotID')['DiseaseID']).max(level=0)
matrix_prot_dis_.columns

len(set(matrix_prot_dis_.columns))
matrix_prot_dis_.shape # --> (8985, 6071)

matrix_prot_dis = pd.DataFrame(matrix_prot_dis_, columns= matrix_prot_dis_.columns, index= not_isolated_proteins)
matrix_prot_dis = matrix_prot_dis.fillna(int(0))
matrix_prot_dis = matrix_prot_dis.astype(int)
matrix_prot_dis.head(2)
matrix_prot_dis.to_csv(os.path.join(wdir, 'mat_protein_dis.txt'), index=False, header=False, sep=" ") 



# DTI  (DRUG - PROTEIN)
matrix_drug_protein_ = pd.get_dummies(DTI.set_index('Drug')['Protein']).max(level=0)
matrix_drug_protein_.columns

len(set(matrix_drug_protein_.columns))
matrix_drug_protein_.shape # -->  (7627, 4884)

matrix_drug_protein = pd.DataFrame(matrix_drug_protein_, columns= not_isolated_proteins, index= not_isolated_drugs)
matrix_drug_protein = matrix_drug_protein.fillna(int(0))
matrix_drug_protein = matrix_drug_protein.astype(int)
matrix_drug_protein.head(2)
matrix_drug_protein.to_csv(os.path.join(wdir, 'mat_drug_protein.txt'), index=False, header=False, sep=" ") 

