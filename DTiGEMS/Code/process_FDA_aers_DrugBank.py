import os
import re
import logging
import argparse
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
from collections import Counter
import xml.etree.ElementTree as ET
import helper_functions_dtigems as hf

######################################## START MAIN #########################################
#############################################################################################
def main():
	# get the parameters from the user
	parser = argparse.ArgumentParser()
	parser.add_argument("dbPath", help="Path to the database",
						type=str)
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
						help="Verbosity (between 1-4 occurrences with more leading to more "
							"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
							"DEBUG=4")
	args = parser.parse_args()
	# define logging level
	log_levels = {
		0: logging.CRITICAL,
		1: logging.ERROR,
		2: logging.WARN,
		3: logging.INFO,
		4: logging.DEBUG,
	}

	# set the logging info
	level= log_levels[args.verbosity]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)

	# sanity check for the DB
	paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
	logging.info(f'\n{paper_cite}\n')
	DB_PATH = args.dbPath
	# DB_PATH =  './../../DB/Data/DrugBank/DrugBank_DTIs.tsv'
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	START_YEAR = 2004
	END_YEAR   = 2021
	hf.dowload_FAERS(START_YEAR,END_YEAR)
	hf.uncompress_FAERS('./../../DB/Data/cross_side_information_DB/FDA')
	all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
	
	####################### DAVIS -- DB annotation specific ###########################
	# create the dicctionary with the drugs up to date for BINDINGDB-DrugBank
	annotation_file = f'./../../DB/Data/cross_side_information_DB/FDA/DrugBank_Drug_FDA_dict{START_YEAR}_{END_YEAR}.txt'
	if os.path.isfile(annotation_file):
		logging.info(f'Found {annotation_file}')
		fda_DB_dict = hf.read_annotation(annotation_file)
	else: 
		all_entries = list(set(hf.get_drugs_DrugBank(DB_PATH)))
		# get the names from drugbank
		tree = ET.parse('./../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
		root = tree.getroot()
		drugbank_dic = {} # dict with generic or product name as key and drugbank_id as value
		for drug_entry in tqdm(root):
			drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
			if drugbank_ID in all_entries:
				name = drug_entry.find('{http://www.drugbank.ca}name').text.upper()
				prod = drug_entry.find('{http://www.drugbank.ca}products')
				prod_names = set([brandname.find('{http://www.drugbank.ca}name').text.upper() for brandname in prod])
				if name:
					drugbank_dic[name] = drugbank_ID
				if len(prod_names) >= 1:
					for prod in prod_names:
						drugbank_dic[prod] = drugbank_ID
		# get the FDA drug names to DrugBank ID
		blacklist = ['TABLETS', 'CAPSULES', 'INJECTABLE', 'INJECTION', 'MG', 'VIAL', 'ML']
		fda_DB_dict = []
		for fda_id, drug_name in tqdm(all_drugs, desc='Crossing IDs with the FDA'):
			if drugbank_dic.get(drug_name, None):
				fda_DB_dict.append((fda_id, drug_name, drugbank_dic.get(drug_name)))
			else:
				drug_name_cleaned = hf.clean_name(blacklist, drug_name)
				if drugbank_dic.get(drug_name_cleaned, None):
					fda_DB_dict.append((fda_id, drug_name_cleaned, drugbank_dic.get(drug_name_cleaned)))
		# write the dictionary
		fda_DB_dict = list(set(fda_DB_dict))
		assert len(set([len(event) for event in fda_DB_dict])) == 1, 'The FDA annotation is not well formatted'
		with open(annotation_file, 'w') as f:
			_ =f.write('fda_incident, drugname, DB_id\n')
			for entry in fda_DB_dict:
				_ = f.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\n')

	############################################################################
	# get the FDA event codes and adverse events
	all_events =  hf.get_events(START_YEAR, END_YEAR)

	# get the events for ONLY the kept drugs
	ids_2_keep = set([drug_entry[0]  for drug_entry in all_drugs])
	all_events = [event for event in tqdm(all_events, desc='Removing drug entries missing in dataset') if event[0] in ids_2_keep]

	# get the keywords from the events
	all_keywords = [element for event in all_events for element in event]
	all_keywords = [element for element in tqdm(all_keywords, desc= 'Extracting IDs') if bool(re.search(r'\D',element))]

	freq = Counter(all_keywords)
	logging.info(
		''' Keywords that appear too frequently (>0.001) or too rarely (the number of reports were <5) were removed.
		remove the reports were <5
		CITE: 
		Takarabe, Masataka, et al. "Drug target prediction using adverse event report systems: 
		a pharmacogenomic approach." Bioinformatics 28.18 (2012): i611-i618.
		''')
	freq_percentage = { drug_id: envents/len(all_keywords) for drug_id, envents in freq.items() if envents > 5}
	freq_percentage = { drug_id: envents for drug_id, envents in freq_percentage.items() if envents < 0.001}
	unique_keywords= set(freq_percentage.keys())
	
	drugs =  hf.get_drugs_DrugBank(DB_PATH)
	drugs = list(set(drugs))
	unique_keywords_vector = sorted(list(unique_keywords))

	aers_freq, aers_bit = hf.get_frequencies(fda_DB_dict, 
											all_events, 
											freq_percentage, 
											drugs, 
											unique_keywords_vector)

	path = hf.check_and_create_folder(db_name)
	aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
	aers_freq_pd.to_csv(os.path.join(path, 'DrugBank_Drug_FDA_aers_freq.tsv'), sep='\t')
	aers_bit_pd =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector, index=drugs)
	aers_bit_pd.to_csv(os.path.join(path, 'DrugBank_Drug_FDA_aers_bit.tsv'), sep='\t')



#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################


