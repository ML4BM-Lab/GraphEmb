import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
from pubchempy import Compound


# DTINet uses DrugBank for Drug Nodes, DTI and Drug-Drug Interaction
# output now in: /home/uveleiro/data/uveleiro/tmp_test_data_DTI

# read Database xml
tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
root = tree.getroot()

##############################################################
################ DRUG NODES - DrugBank #######################

# output files: drugs.txt & drug_dic_map.txt
# with columns: DB_ID/KEGG_ID/SMILES/name

# DTINet does not need information about KEGG compound/drug number
# return a tuple that contanis only drugs ID & name: (drugbank_ID, name).
# this can be modified as in "process_DrugBank_getKEGG.py"

#### Here we return the initial files for drugs.
# Returning only DB idenfitier and name
db_drug_nodes = []
for drug_entry in tqdm(root):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	name = drug_entry.find('{http://www.drugbank.ca}name').text # incluir aqui name sin check
	drug = (drugbank_ID, name) # crea tupla drug ID number & name
	db_drug_nodes.append(drug)


# Drug nodes (drugs.txt). All nodes in DB (w\o filter) 
with open('/home/uveleiro/data/uveleiro/tmp_test_data_DTI/drugs.txt', 'w') as f:
	#_ = f.write('#DrugBankID\n')
	for item in db_drug_nodes:
		_ = f.write("%s\n" % item[0])

# Write dictionary (drugbank_ID and name)
with open('/home/uveleiro/data/uveleiro/tmp_test_data_DTI/drug_dic_map.txt', 'w') as f:
	for item in db_drug_nodes:
		_ = f.write("%s:%s\n" % (item[0],item[1]))



### Get SMILES 
# (needed for similarity measures with fingerprints!)

def get_compound_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

db_drug_smiles  = []
for drug_entry in tqdm(root):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	smiles = 'None'
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
				break
	if smiles == 'None':
		for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
			for ids in exids:
				if(ids.text == 'PubChem Compound'): 
					pubchem_id = exids[1].text
					smiles = get_compound_pubchem(pubchem_id)
					break
	drug = (drugbank_ID, smiles)
	db_drug_smiles.append(drug)

with open('/home/uveleiro/data/uveleiro/tmp_test_data_DTI/drugs_w_SMILES.tsv', 'w') as f:
	#_ = f.write('# DrugBank ID\tProtein list\n')
	for item in db_drug_smiles:
		_ = f.write("%s\t%s\n" % item)



###########################################################################
######################## DTIs - DrugBank ##################################

def get_targets_list(drug_entry):
	'''
	This function gets a drug entry and returns a tuple that contains
	the drugbankID and a list of interacting targets
	uses list() instead of .getchildren(), the second will be removed in future versions
	'''
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
with open('/home/uveleiro/data/uveleiro/tmp_test_data_DTI/relation_drug_prot_list.tsv', 'w') as f:
	#_ = f.write('# DrugBank ID\tProtein list\n')
	for item in db_rel_drug_prot:
		_ = f.write("%s\t%s\n" % item)

###################################################################

#### GET COORDINATES - FUNCTION
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

## save file
with open('/home/uveleiro/data/uveleiro/tmp_test_data_DTI/DTI_coordinates.tsv', 'w') as f:
	for item in coordinate_list:
		_ = f.write("%s\t%s\n" % item)


## aquí faltaria crear la matriz
####

####
# usar funcion comun

###########################################################
######## Drug-Drug Interactions - DrugBank ################

def get_drug_drug_coordinates(drug_entry, drug_drug_coodinates):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	all_interactions =  drug_entry.findall('.//{http://www.drugbank.ca}drug-interactions')[0]
	for inter in all_interactions: # aqui iteramos en targets
		coordinate = (drugbank_ID, list(inter)[0].text)	
		#print(coordinate)
		drug_drug_coodinates.append(coordinate)
	return drug_drug_coodinates

drug_drug_coodinates = []
for drug_entry in tqdm(root):
	drug_drug_coodinates = get_drug_drug_coordinates(drug_entry, drug_drug_coodinates)

with open('/home/uveleiro/data/uveleiro/tmp_test_data_DTI/drug-drug_coordinates.tsv', 'w') as f:
	for item in drug_drug_coodinates:
		_ = f.write("%s\t%s\n" % item)