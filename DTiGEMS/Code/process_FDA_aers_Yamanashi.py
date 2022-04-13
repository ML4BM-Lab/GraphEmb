
import os, sys
import re
import logging
import pandas as pd
import zipfile
import urllib.request
import requests as rq
from tqdm import tqdm
from collections import Counter
import multiprocessing as mp
from functools import reduce
import argparse
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

	# parse the ascii ones not the xmls

	# set the logging info
	level= log_levels[args.verbosity]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)

	# sanity check for the DB
	paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
	logging.info(f'\n\n{paper_cite}\n')
	DB_PATH = args.dbPath

	# DB_PATH =  './../../DB/Data/Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt'
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
		
	START_YEAR = 2004
	END_YEAR = 2021

	hf.dowload_FAERS(START_YEAR,END_YEAR)
	hf.uncompress_FAERS('./../../DB/Data/cross_side_information_DB/FDA')

	all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
	
	####################### Yamanashi -- KEGG annotation specific ###########################
	# create the dicctionary with the drugs up to date for KEGG
	annotation_file = f'./../../DB/Data/cross_side_information_DB/FDA//fda_kegg_dict{START_YEAR}_{END_YEAR}.txt'
	if os.path.isfile(annotation_file):
		logging.info(f'Found {annotation_file}, reading it')
		fda_kegg_dict = hf.read_annotation(annotation_file)
		assert len(set([len(event) for event in fda_kegg_dict])) == 1, 'The FDA annotation is not well formatted'
	else:
		KEGG_NAME_DB = re.compile(r'(^[^\(\d{2}\)]*) ')
		PATTERN_DB_KEGG = re.compile(r"(?<=DR:)[\w]+")
		kegg_drug_names = rq.get('http://rest.kegg.jp/list/drug/').text.split('\n')
		kegg_drug_names = [entry.upper() for entry in kegg_drug_names if entry]
		kegg_dic = {}
		for entry in tqdm(kegg_drug_names):
			kegg_id, kegg_desc = entry.split('\t')
			parsed_desc = [KEGG_NAME_DB.findall(pseudo.strip()) for pseudo in kegg_desc.split(';')]
			parsed_desc = [item for sublist in parsed_desc for item in sublist]
			for name in parsed_desc:
				kegg_dic[name] = kegg_id.strip()
		#
		blacklist = ['TABLETS', 'CAPSULES', 'INJECTABLE', 'INJECTION', 'MG', 'VIAL', 'ML']
		fda_kegg_dict = []
		for fda_id, drug_name in tqdm(all_drugs):
			kegg_id = kegg_dic.get(drug_name, None)
			if kegg_id:
				fda_kegg_dict.append((fda_id.strip(), drug_name.strip(), kegg_id.strip()))
			else:
				drug_name_cleaned = re.sub(r'\([^)]*\)', '', drug_name)
				drug_name_cleaned = re.sub("[(,)]","", drug_name_cleaned)
				drug_name_cleaned = re.sub("[\d]+","", drug_name_cleaned)
				drug_name_cleaned = [word for word in drug_name_cleaned.split() if word not in blacklist]
				drug_name_cleaned = ' '.join(drug_name_cleaned).strip()
				kegg_id = kegg_dic.get(drug_name_cleaned, None)
				if kegg_id:
					kegg_id = PATTERN_DB_KEGG.findall(kegg_id)[0].strip()
					fda_kegg_dict.append((fda_id.strip(), drug_name_cleaned.strip(), kegg_id))
		fda_kegg_dict = list(set(fda_kegg_dict))
		assert len(set([len(event) for event in fda_kegg_dict])) == 1, 'The FDA annotation is not well formatted'
		with open(annotation_file, 'w') as f:
			_ =f.write('fda_incident, drug_name, kegg_id\n')
			for entry in tqdm(fda_kegg_dict):
				_ = f.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\n')
	############################################################################
	all_events =  hf.get_events(START_YEAR, END_YEAR)

	# get the events for ONLY for the kept drugs (PS and SS)
	# beware with the Ids comming from later events, have two different number cases
	ids_2_keep = set([drug_entry[0]  for drug_entry in all_drugs])
	all_events = [event for event in tqdm(all_events) if event[0] in ids_2_keep]

	# get the keywords from the events
	all_keywords = [element for event in all_events for element in event]
	all_keywords = [element for element in tqdm(all_keywords) if bool(re.search(r'\D',element))]

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
	unique_keywords = set(freq_percentage.keys())

	drugs = hf.read_and_extract_drugs(DB_PATH)
	drugs = list(set(drugs))

	unique_keywords_vector = sorted(list(unique_keywords))

	#IdDrug dict
	fda_kegg_dict_drug_id = {}
	for event_code, _, kegg_id in tqdm(fda_kegg_dict):
		if kegg_id:
			if kegg_id in fda_kegg_dict_drug_id:
				fda_kegg_dict_drug_id[kegg_id].append(event_code)
			else:
				fda_kegg_dict_drug_id[kegg_id] = [event_code]

		#IDEvents AERS dict
	id_events_aers = {}
	for idn, aersname in tqdm(all_events):
		if idn in id_events_aers:
			id_events_aers[idn].append(aersname)
		else:
			id_events_aers[idn] = [aersname]

	aers_freq = []
	for drug in tqdm(drugs):
		aers_freq.append(hf.get_freqs(drug, 
									unique_keywords_vector, 
									fda_kegg_dict_drug_id, 
									freq_percentage, id_events_aers))

	aers_bit = []
	for freq in aers_freq:
		aers_bit.append([1 if entry > 0 else 0 for entry in freq])

	aers_bit =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector,index=drugs)
	aers_freq = pd.DataFrame(aers_freq, columns=unique_keywords_vector,index=drugs)
	path = hf.check_and_create_folder(db_name)
	subname = hf.get_yamanashi_subDB(db_name)
	aers_freq.to_csv(os.path.join(path, f'{subname}_Drug_FDA_aers_freq.tsv'), sep='\t')
	aers_bit.to_csv(os.path.join(path,  f'{subname}_Drug_FDA_aers_bit.tsv'), sep='\t')


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################