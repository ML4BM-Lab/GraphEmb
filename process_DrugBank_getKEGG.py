from pickle import NONE
import xml.etree.ElementTree as ET
from tqdm import tqdm

tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
root = tree.getroot()

# WARNING: this only retrieves the KEGG drug IDs, not the Compound IDs as far as
# the SIMCOMP score only uses drugs not compounds...


def check_drug(entry):
	try:
		kegg_index = [ _.text for _ in entry.findall('.//{http://www.drugbank.ca}resource')].index('KEGG Drug')
		kegg_ID = entry.findall('.//{http://www.drugbank.ca}identifier')[kegg_index].text
		return kegg_ID
	except ValueError:
		return None

def check_compound(entry):
	try:
		kegg_index = [ _.text for _ in entry.findall('.//{http://www.drugbank.ca}resource')].index('KEGG Compound')
		kegg_ID = entry.findall('.//{http://www.drugbank.ca}identifier')[kegg_index].text
		return kegg_ID
	except ValueError:
		None

def check_compound_smile_name(entry):
	try:
		kegg_index = [ _.text for _ in entry.findall('.//{http://www.drugbank.ca}resource')].index('KEGG Compound')
		kegg_ID = entry.findall('.//{http://www.drugbank.ca}identifier')[kegg_index].text
		smile = None
		name = None
		return((drugbank_ID, kegg_ID, smile, name))
	except ValueError:
		None


db_kegg_drugs     = []
db_kegg_compounds = []
for drug_entry in tqdm(root):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	drug = (drugbank_ID, check_drug(drug_entry))
	if drug[1]:
		db_kegg_drugs.append(drug)
	compound = (drugbank_ID, check_compound(drug_entry))
	if compound[1]:
		db_kegg_compounds.append(compound)

with open('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/DB_2_KEGG_drugs.tsv', 'w') as f:
	_ = f.write('# DrugBank ID\tKEGG Drug ID\n')
	for item in db_kegg_drugs:
		_ = f.write("%s\t%s\n" % item)

with open('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/DB_2_KEGG_drugsCompounds.tsv', 'w') as f:
	_ = f.write('# DrugBank ID\tKEGG Compound ID\n')
	for item in db_kegg_compounds:
		_ = f.write("%s\t%s\n" % item)




for index in range(0, len(drug_entry)):
		print(f'{index=}, {drug_entry[index].tag=}, {drug_entry[index].text=}')


		index.findall('.//{http://www.drugbank.ca}drug-interactions')


results = []
for drug_entry in tqdm(root):
	for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
		for prop in props:
			if(prop.text == 'SMILES'):
				results.append((props[0].text, props[1].text))
