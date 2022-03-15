import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
from pubchempy import Compound
import json
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
from itertools import repeat

# change output to (relative path): 
output_path = '../Data/DrugBank'

# read xml
tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
root = tree.getroot()

##############################################################
################ DRUG NODES - DrugBank #######################
# output files: 
#     - all_drugs.txt
#     - all_drug_dic_map.txt
#     - dic_drugnames_DBID.json

db_drug_nodes = []
drug_IDs = []
drug_names = []
for drug_entry in tqdm(root):
    drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
    name = drug_entry.find('{http://www.drugbank.ca}name').text # incluir aqui name sin check
    drug_names.append(name)
    drug_IDs.append(drugbank_ID)

# Drug nodes (all_drugs.txt). All nodes in DB (w\o filter) 
with open(os.path.join(output_path, 'all_drugs.txt'), 'w') as f:
	#_ = f.write('#DrugBankID\n')
	for item in drug_IDs:
		_ = f.write("%s\n" % item)

# Write dictionary (drugbank_ID and name)
with open(os.path.join(output_path,'all_drug_dic_map.txt'), 'w') as f:
	for i in range(len(drug_IDs)):
		_ = f.write("%s:%s\n" % (drug_IDs[i],drug_names[i]))

#  ------- uncomment to save it as JSON too ---------
# dic_drugnames_DBID = dict(zip(drug_IDs, drug_names)) 
# save this dictionary in cross_side..._DB
# file_name_json = 'dic_drugnames_DBID.json'
# with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
#     json.dump(dic_drugnames_DBID, f, ensure_ascii=False, indent=4)


###########################################################################
######################## DTIs - DrugBank ##################################
'''
# written but not used for DTINet
def get_targets_list(drug_entry):
	#This function gets a drug entry and returns a tuple that contains
	#the drugbankID and a list of interacting targets
	#uses list() instead of .getchildren(), the second will be removed in future versions
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	target_list = []
	for tgt in drug_entry.findall('.//{http://www.drugbank.ca}targets')[0]: # iterate in targets
		for i in tgt.findall('.//{http://www.drugbank.ca}polypeptide'): # for searching external-identifiers
			ext = i.findall('.//{http://www.drugbank.ca}external-identifiers')[0]
			for id in list(ext):
				if list(id)[0].text == "UniProtKB": 
					uniprot_ID = list(id)[1].text 
					target_list.append(uniprot_ID)
	db_line = (drugbank_ID, target_list)
	return db_line

## execute
db_rel_drug_prot = []
for drug_entry in tqdm(root):
	line_drug_target = get_targets_list(drug_entry)
	db_rel_drug_prot.append(line_drug_target)

## save
with open(os.path.join(output_path ,'relation_drug_prot_list.tsv'), 'w') as f:
	#_ = f.write('# DrugBank ID\tProtein list\n')
	for item in db_rel_drug_prot:
		_ = f.write("%s\t%s\n" % item)
'''


#### GET COORDINATES of DTI INTERACTIONS ##########
def get_DTI_coordinates(drug_entry, coordinate_list):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	target = 'None'
	for tgt in  drug_entry.findall('.//{http://www.drugbank.ca}targets')[0]: # iterate in targets
		for i in tgt.findall('.//{http://www.drugbank.ca}polypeptide'): # for searching external-identifiers
			ext = i.findall('.//{http://www.drugbank.ca}external-identifiers')[0]
			for id in list(ext):
				if list(id)[0].text == "UniProtKB": # busca donde esta uniprot
					target = list(id)[1].text
					if target != 'None':
						coordinate = (drugbank_ID, target)
						coordinate_list.append(coordinate)
					else:
						continue 
	return coordinate_list

## execute the function
coordinate_list = []
for drug_entry in tqdm(root):
	coordinate_list = get_DTI_coordinates(drug_entry, coordinate_list)
	#db_rel_drug_prot.append(line_drug_target)

## save coordinates to a tsv
with open(os.path.join(output_path, 'coordinates_DTI.tsv'), 'w') as f:
	for item in coordinate_list:
		_ = f.write("%s\t%s\n" % item)

## GET MATRIX 
#corrdinate_list_unique = list(set(coordinate_list)) # unique interactiosn; faster no duplicity (changes order but not so important here)
# using the list with repeated coordinates, because it mantains the order, then we get Drug_ID starting with DB00001
df_t = pd.DataFrame(coordinate_list, columns=['Drug_ID', 'Protein_ID'])
matrix_drug_protein = pd.get_dummies(df_t.set_index('Drug_ID')['Protein_ID']).max(level=0) #.reset_index()
matrix_drug_protein

# index false & header false for final matrix
matrix_drug_protein.to_csv(os.path.join(output_path, 'mat_drug_protein.tsv'), index=True, header=True, sep='\t') 

# for transporse matrix (needed for all files)
# we need protein_drug.txt & mat_drug_protein_remove_homo.txt; but get with the "processed matrix"
# with matrix_protein_drug = matrix_drug_protein.transpose() 


#############################################################################
###################### Drug-Drug Interactions - DrugBank ####################

def get_drug_drug_coordinates(drug_entry, drug_drug_coodinates):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	all_interactions =  drug_entry.findall('.//{http://www.drugbank.ca}drug-interactions')[0]
	for inter in all_interactions: # aqui iteramos en targets
		coordinate = (drugbank_ID, list(inter)[0].text)	
		drug_drug_coodinates.append(coordinate)
	return drug_drug_coodinates

drug_drug_coodinates = []
for drug_entry in tqdm(root):
	drug_drug_coodinates = get_drug_drug_coordinates(drug_entry, drug_drug_coodinates)


#coordinates df
df_d = pd.DataFrame(drug_drug_coodinates, columns=['Drug_ID_A', 'Drug_ID_B']) 
df_d = df_d.drop_duplicates()
df_d.to_csv(os.path.join(output_path, 'coordinates_drug_drug.tsv'), header=False,index=False ,sep="\t")

'''
#matrix
matrix_drug_drug = pd.get_dummies(df_d.set_index('Drug_ID_A')['Drug_ID_B']).max(level=0) #.reset_index()# no esta ordenado esto
matrix_drug_drug

# index false & header false for final matrix
matrix_drug_drug.to_csv(os.path.join(output_path, 'mat_drug_drug.tsv'), index=True, header=True, sep='\t') 
'''