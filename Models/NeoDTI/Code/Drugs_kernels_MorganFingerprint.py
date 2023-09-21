from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import os
from tqdm import tqdm
import numpy as np
import sys
import pandas as pd
from itertools import repeat
import json


def get_drug_similarity_matrix(file_path_simdrug_mat, list_of_drug_nodes, dict_drugid_smiles):
	list_drugs = list_of_drug_nodes
	list_smiles = [dict_drugid_smiles[i] for i in list_drugs] # 
	all_Sim_Tani = []
	dic = {}
	for drug in tqdm(range(len(list_drugs)), desc='Retrieving Pairwise Tanimoto for drugs', position=0, leave=True): 
		id_drug_1 = list_drugs[drug]
		smiles_drug_1 = list_smiles[list_drugs.index(id_drug_1)]
		sim_for_drug = []
		tmp = []
		tmp.extend(repeat(smiles_drug_1, len(list_drugs))) 
		for j in range(len(tmp)):
			result = get_pairwise_dice(tmp[j], list_smiles[j], dic)
			sim_for_drug.append(result)
		all_Sim_Tani.append(sim_for_drug)
	#
	df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= list_drugs, index = list_drugs) # add columns, & index
	# save pickle
	#df_all_Sim_Tani.to_pickle(file_path_sim_pickle) # add here column %& index next time
	# save csv for model
	df_all_Sim_Tani.to_csv(file_path_simdrug_mat, sep='\t', header=False, index=False) # add here column %& index next time


def get_pairwise_dice(smiles1,smiles2, dic): #dic == smile2fp
	try:
		for smile in [smiles1, smiles2]:
			if not smile in dic:
				mol1 = Chem.MolFromSmiles(str(smile))
				#fp1  = Chem.RDKFingerprint(mol1)
				fp1 = AllChem.GetMorganFingerprint(mol1,2) #RADIUS 2
				dic[smile] = fp1
		#tani = DataStructs.FingerprintSimilarity(dic[smiles1],dic[smiles2]) #pairwise similarity
		tani = DataStructs.DiceSimilarity(dic[smiles1],dic[smiles2])
		return tani
	except:
		return None

try:
	subdataset = sys.argv[2]
except:
	subdataset = ""
dataset = sys.argv[1]
wdir = f'NeoDTI/Data/{dataset}/{subdataset}'
file_path_simdrug_mat = os.path.join(wdir, 'Similarity_Matrix_Drugs.txt')
list_of_drug_nodes = np.loadtxt(os.path.join(wdir,'drug.txt'), dtype='str').tolist()
# save json 
with open(os.path.join(wdir,'dic_smiles.json'), 'r') as f:
	dict_drugid_smiles = json.load(f)

if ((not os.path.exists(file_path_simdrug_mat)) ):
	# call function: 
	get_drug_similarity_matrix(file_path_simdrug_mat, list_of_drug_nodes, dict_drugid_smiles)
	#pass
