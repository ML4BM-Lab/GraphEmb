import os
import pandas as pd
from tqdm import tqdm
from pubchempy import Compound
from rdkit import Chem
import xml.etree.ElementTree as ET

def get_compound_pubchem(drug):
	try:
		return Compound.from_cid(drug).isomeric_smiles
	except:
		return None

def get_smiles(drug_entry):
	'''
	Get list of drugs with smiles from DrugBank xml, 
	for SMILES that are not in DrugBank retrieves them from PubChem
	Then, check if we can create a fingerprint (if not, do not include)
	'''
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	smiles = None
	fp = None
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props: 
			if(prop.text == 'SMILES'):
				smiles = props[1].text
				#print("Hi" + smiles)
				#Chem.MolFromSmiles(str(smiles))
				break
	if not smiles:
		for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
			for ids in exids:
				if(ids.text == 'PubChem Substance'): 
					pubchem_id = exids[1].text
					smiles = get_compound_pubchem(pubchem_id)
					#print(smiles)
					break
	if not smiles:
		return(drugbank_ID, None)
	else:
		return(drugbank_ID, smiles)

def get_drug_smiles_drugbank(output_folder):
	tree = ET.parse(os.getcwd()+ '/../../DB/Data/DrugBank/full_database.xml')
	root = tree.getroot()
	list_drugs  = []
	list_smiles = []
	for drug_entry in tqdm(root):
		drug_id, smiles = get_smiles(drug_entry)
		if drug_id and smiles:
			list_drugs.append(drug_id)
			list_smiles.append(smiles) 
	assert len(list_drugs) == len(list_smiles), 'The length of the Drug IDs does not match the number of SMILES'
	df = pd.DataFrame()
	df['DrugBank_ID'] = list_drugs
	df['SMILES'] = list_smiles

	output_path = os.getcwd() + '/../Data/'+output_folder+'/drugs_smiles.txt'
	df.to_csv(output_path, header=None, index = None, sep = ' ')