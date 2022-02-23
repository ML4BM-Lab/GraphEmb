import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd
import numpy as np
import os
from pubchempy import Compound


# este script es el bueno!!! no borrar!!!!
##### falta drug-drug!!

# all output now in uveleiro/Test

# DTINet uses DrugBank for Drug Nodes and DTI 

### DRUG NODES - DrugBank 
# output files: info_drugs.tsv & info_compounds.tsv
# with columns: DB_ID/KEGG_ID/SMILES/name

# first part is similar to "process_DrugBank_getKEGG.py"

tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
root = tree.getroot()


###############################################################
###############################################################

# DTINet does't need information about KEGG compound/drug number
# return a tuple of drugs (drugbank_ID, name).

#### Here we return the initial files for drugs.
# Returning only DB idenfitier and name
db_drug_nodes = []
for drug_entry in tqdm(root):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	name = drug_entry.find('{http://www.drugbank.ca}name').text # incluir aqui name sin check
	drug = (drugbank_ID, name) # crea tupla drug ID number & name
	db_drug_nodes.append(drug)


## writting files
# write the initial file that contain drug nodes 
with open('/home/uveleiro/data/uveleiro/Test/test_output/drug_initial_file.txt', 'w') as f:
	#_ = f.write('#DrugBankID\n')
	for item in db_drug_nodes:
		_ = f.write("%s\n" % item[0])

# write dictionary file
with open('/home/uveleiro/data/uveleiro/Test/test_output/drug_dic_map.txt', 'w') as f:
	for item in db_drug_nodes:
		_ = f.write("%s:%s\n" % (item[0],item[1]))


# this part is similar to what we have above, but returning the SMILES
db_drug_smiles = []
for drug_entry in tqdm(root):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	# ---- bucle para sacar info de los smiles; poner más bonito en un futuro si se puede -------
	smiles ='None'
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
	# ---- acaba aqui------
	drug = (drugbank_ID, smiles) # crea tupla drug ID number & name
	db_drug_smiles.append(drug)

for item in db_drug_smiles[:10]:
	print(item, end='\n')

df_drugs = pd.DataFrame(db_drug_smiles)
df_drugs.columns = ['ID', 'SMILES']
df_drugs.SMILES[df_drugs.SMILES == "None"].count() 
df_drugs.shape

with open('/home/uveleiro/data/uveleiro/Test/test_output/drugs_w_SMILES.tsv', 'w') as f:
	#_ = f.write('# DrugBank ID\tProtein list\n')
	for item in db_drug_smiles:
		_ = f.write("%s\t%s\n" % item)

# hay bastante diferencia entre el numero de nodos que tienen ellos 
# y lo que da aquí (aun quitando los smiles que faltan)

### no todos tienen smiles asi que hay que arreglar esta funcion
#### volver sobre esto porque faltan smiles de algunos compuestos
#### => usar pubchem id y sacarlo de ahí



def get_drug_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

db_drug_smiles  = []
for drug_entry in tqdm(root):
	#drug_entry = root[0]
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	#drugbank_ID
	# ---- bucle para sacar info de los smiles; poner más bonito en un futuro si se puede -------
	smiles = 'None'
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
				break
		if smiles == 'None':
			for props in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
				for prop in props:
					if(prop.text == 'PubChem Substance'): 
						#print(props[1].text)
						pubchem_id = props[1].text
						smiles = get_drug_pubchem(pubchem_id)
					break # romper una vez que encuentre el pubchem
	drug = (drugbank_ID, smiles) # crea tupla drug ID number & name
	db_drug_smiles.append(drug)

# export results to a tsv file
# realmente lo mejor seria hacer el write dentro del bucle para que no pese tanto en memoria!
# cambiar! pero esto funciona!
with open('/home/uveleiro/data/uveleiro/Test/test_output/drugs_w_SMILES_complete.tsv', 'w') as f:
	#_ = f.write('# DrugBank ID\tProtein list\n')
	for item in db_drug_smiles:
		_ = f.write("%s\t%s\n" % item)

###########################################################################
##########################################################################


# funcion que devuelve una tabla
def get_targets_list(drug_entry):
	'''
	This function gets a drug entry and returns a tuple that contains
	the drugbankID and a list of interacting targets
	uses list() instead of .getchildren(), the second will be removed in future versions
	'''
	drug_name = drug_entry.find('{http://www.drugbank.ca}name').text
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	target_list = []
	all_targets =  drug_entry.findall('.//{http://www.drugbank.ca}targets')[0]
	for tgt in all_targets: # aqui iteramos en targets
		tgt.find('{http://www.drugbank.ca}name').text  # print el nombre para comprobar
		# entrar en polypeptide
		plp = tgt.findall('.//{http://www.drugbank.ca}polypeptide')
		for i in plp:
			# encontrar los external identifiers
			ext = i.findall('.//{http://www.drugbank.ca}external-identifiers')[0]
			for id in list(ext):#.getchildren():
				if list(id)[0].text == "UniProtKB": # busca donde esta uniprot
					uniprot_ID = list(id)[1].text # selecciona uniprot
					target_list.append(uniprot_ID)
	db_line = (drugbank_ID, target_list)
	return db_line

#### execute
db_rel_drug_prot = []
for drug_entry in tqdm(root):
	line_drug_target = get_targets_list(drug_entry)
	db_rel_drug_prot.append(line_drug_target)


## escribir tabla lista
with open('/home/uveleiro/data/uveleiro/Test/test_output/relation_drug_prot_list.tsv', 'w') as f:
	#_ = f.write('# DrugBank ID\tProtein list\n')
	for item in db_rel_drug_prot:
		_ = f.write("%s\t%s\n" % item)

#################
#################

# Pero realmente queremos una lista de coordenadas
# tuplas tal que (DB_ID, UniprotID)


##### GET COORDINATES - FUNCTION
def get_DTI_coordinates(drug_entry, coordinate_list):
	# drug_name = drug_entry.find('{http://www.drugbank.ca}name').text
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	target = 'None'
	#coordinate_list = []
	all_targets =  drug_entry.findall('.//{http://www.drugbank.ca}targets')[0]
	for tgt in all_targets: # aqui iteramos en targets
		plp = tgt.findall('.//{http://www.drugbank.ca}polypeptide')
		for i in plp:
			ext = i.findall('.//{http://www.drugbank.ca}external-identifiers')[0]
			for id in list(ext):
				if list(id)[0].text == "UniProtKB": # busca donde esta uniprot
					target = list(id)[1].text
					if target != 'None':
						coordinate = (drugbank_ID, target)
						#print(coordinate)
						coordinate_list.append(coordinate)
					else:
						continue 
	return coordinate_list

## execute the function
coordinate_list = []
for drug_entry in tqdm(root):
	coordinate_list = get_DTI_coordinates(drug_entry, coordinate_list)
	#db_rel_drug_prot.append(line_drug_target)

## write a coordinate file
## escribir tabla lista
with open('/home/uveleiro/data/uveleiro/Test/test_output/DTI_coordinates.csv', 'w') as f:
	for item in coordinate_list:
		_ = f.write("%s\t%s\n" % item)

# test
df_coordinates = pd.DataFrame(coordinate_list)
df_coordinates.columns = ['DRUG_ID', 'TARGET_ID']
df_coordinates

## aquí faltaria crear la matriz


###########################################################
###########################################################
# Drug-Drug Interactions

# Parecido a lo anterior pero tiene que estar en otro sitio en vez de en targets

def get_drug_drug_coordinates(drug_entry, drug_drug_coodinates):
	#drug_entry = root[55]
	#drug_name = drug_entry.find('{http://www.drugbank.ca}name').text
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	#drugbank_ID
	all_interactions =  drug_entry.findall('.//{http://www.drugbank.ca}drug-interactions')[0]
	for inter in all_interactions: # aqui iteramos en targets
		coordinate = (drugbank_ID, list(inter)[0].text)	
		#print(coordinate)
		drug_drug_coodinates.append(coordinate)
	return drug_drug_coodinates

drug_drug_coodinates = []
for drug_entry in tqdm(root):
	drug_drug_coodinates = get_drug_drug_coordinates(drug_entry, drug_drug_coodinates)

with open('/home/uveleiro/data/uveleiro/Test/test_output/drug-drug_coordinates.csv', 'w') as f:
	for item in drug_drug_coodinates:
		_ = f.write("%s\t%s\n" % item)
