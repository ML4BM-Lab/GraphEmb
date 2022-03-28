
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
import xml.etree.ElementTree as ET


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['bindingDB', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR']
	for db in DB_NAMES:
		if re.search(db, path):
			logging.info(f'Database: {db}')
			if db in ['E', 'GPCR', 'IC', 'NR']:
				db = os.path.join('Yamanashi_et_al_GoldStandard', db)
				return db
			else:
				return db
	logging.error(f'Database: {db} not found')
	sys.exit('Please provide a valid database')

def read_and_extract_drugs(path):
	drugs = []
	with open(path, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				drug, _ = line.split('\t')
				drugs.append(drug)
	return drugs

def get_drugs_bindingDB(path):
	with open(path, 'r') as f:
		_ = next(f)
		drugs = f.readlines()
	drugs = [entry.strip().split('\t')[6] for entry in drugs]
	return drugs

def dowload_FAERS( start_year = 2004, end_year = 2021):
	"""
	Downloads the FAERS data from the FDA website
	"""
	logging.info(f'Downloading FAERS data from {start_year} to {end_year}')
	# it was a change on 2012Q4 from aers to faers
	for year in range( start_year, end_year+1):
		for quarter in range(1, 5):
			dets_file = f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}.zip'
			if(not os.path.exists(dets_file)):
				logging.info(f'Downloading FAERS data from {year}-Q{quarter}')
				url = f'https://fis.fda.gov/content/Exports/faers_ascii_{year}Q{quarter}.zip'
				logging.debug(f'Downloading {url}')
				try :
					urllib.request.urlretrieve(url, dets_file)
				except urllib.error.HTTPError:
					logging.warning(f'falling back to aers instead of faers for {year}-Q{quarter}')
					url = f'https://fis.fda.gov/content/Exports/aers_ascii_{year}Q{quarter}.zip'
					urllib.request.urlretrieve(url, dets_file)
			else:
				logging.debug(f'Already existing file for {year}-Q{quarter}. Skipping')

def uncompress_FAERS(path_2_folder):
	"""
	Uncompress the FAERS data
	"""
	zipfiles = filter(lambda file: file.endswith('zip') ,os.listdir(path_2_folder))
	for file in zipfiles:
		dest = os.path.join(path_2_folder, file.split('.')[0])
		if os.path.exists(dest):
			logging.debug(f'{dest} already exists. Skipping')
			continue
		os.mkdir(dest)
		logging.info(f'Uncompressing {file}')
		with zipfile.ZipFile(os.path.join(path_2_folder , file), 'r') as zip_ref:
			zip_ref.extractall(dest)

def get_events( start_year = 2004, end_year = 2021):
	all_events = []
	for year in tqdm(range(start_year,end_year+1)):
		for quarter in range(1, 5):
			logging.debug(f'Processing FAERS data from {year}-Q{quarter}')
			# someone screwed up the naming of the files on the FDA side....
			possible_files = [
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/REAC{str(year)[-2:]}Q{quarter}.txt',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/REAC{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/REAC{str(year)[-2:]}Q{quarter}.txt',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/REAC{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.TXT',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.txt'
				]
			fl = [real_file for real_file in possible_files if os.path.isfile(real_file)]
			if fl:
				fl = fl[0]
			else:
				logging.warning(f'No file found for {year}-Q{quarter}')
				continue
			logging.debug(f'Reading {fl}')
			with open(fl, 'r') as f:
				_ = next(f)
				events = f.readlines()
			events_cleaned = []
			# later files from the fda, have 3 fields instead of 2
			for event in events:
				event = event.strip().split('$')
				if len(event) == 3:
					events_cleaned.append([event[0], event[1].upper()])
				else:
					events_cleaned.append([event[0], event[2].upper()])
			all_events.extend(events_cleaned)
			logging.debug(f'{len(all_events)} total events')
	return all_events

def get_drugs( start_year = 2004, end_year = 2021):
	all_drugs = []
	for year in tqdm(range(start_year,end_year+1)):
		for quarter in range(1, 5):
			logging.debug(f'Processing FAERS data from {year}-Q{quarter}')
			# someone screwed up the naming of the files on the FDA side....
			possible_files = [
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/DRUG{str(year)[-2:]}Q{quarter}.txt',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/DRUG{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/DRUG{str(year)[-2:]}Q{quarter}.txt',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/DRUG{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/drug{str(year)[-2:]}q{quarter}.TXT',
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/drug{str(year)[-2:]}q{quarter}.txt'
			]
			fl = [real_file for real_file in possible_files if os.path.isfile(real_file)]
			if fl:
				fl = fl[0]
			else:
				logging.warning(f'No file found for {year}-Q{quarter}')
				continue
			with open(fl, 'r', encoding = "ISO-8859-1") as f:
				_ = next(f)
				drugs = f.readlines()
			all_drugs.extend(drugs)
	return all_drugs

def read_bindingDB_annotation(annotation_file):
	with open(annotation_file, 'r') as f:
		_ = next(f)
		dictionary = f.readlines()
	dictionary = [entry.strip().split('\t') for entry in dictionary]
	return dictionary

def parse_drugs(entry):
	drugs_to_keep = ['PS', 'SS']
	entry = entry.strip().split('$')
	try:
		if entry[2] in drugs_to_keep :
			return (entry[0], entry[3])
		elif entry[3] in drugs_to_keep:
			return  (entry[0], entry[4])
	except IndexError:
		logging.warning(f'{entry} is not a valid drug entry')

def check_and_create_folder(db_name):
	folder =os.path.join('./../Data', db_name)
	if not os.path.exists(folder):
		os.makedirs(folder)
	return folder

def get_freqs(drug, keywords, fda_dict, keyword_freq_percent, id_aers_dict):
	drug_kw_vector = [0]*len(keywords)
	all_cases = fda_dict.get(drug)
	if not all_cases:
		return(drug_kw_vector)
	all_drug_events = [id_aers_dict.get(case,[0]) for case in set(all_cases)]
	all_drug_events = [x for x in all_drug_events if x!=0]
	all_drug_events = [item for sublist in all_drug_events for item in sublist]
	keywords_set = {x:i for i,x in enumerate(keywords)}
	for event in all_drug_events:
		if event in keywords_set:
			drug_kw_vector[keywords_set.get(event)] = keyword_freq_percent.get(event)
	return drug_kw_vector

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
	# DB_PATH =  '/home/margaret/data/jfuente/DTI/Data/BindingDB/BindingDB_All_2021m11/BindingDB_Kd_filtered_columns_subset_HomoSapiens_noduplicities.tsv'
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = get_DB_name(DB_PATH)
		
	START_YEAR = 2004
	END_YEAR = 2021

	dowload_FAERS(START_YEAR,END_YEAR)
	uncompress_FAERS('./../../../Data/cross_side_information_DB/FDA')

	all_drugs = get_drugs(START_YEAR, END_YEAR)
	# unique_drugs = len(set([entry.split('$')[3] for entry in all_drugs]))
	# logging.debug(f'Found {unique_drugs} unique drugs Vs the 291.997 claimed on the paper. ACCEPTABLE')
	logging.info('''
		In this study, we used the drugs (compound names or product names) 
		labeled as PS or SS and searched the equivalent drugs registered in the KEGG database.
		CITE:
		Takarabe, Masataka, et al. "Drug target prediction using adverse event report systems: 
		a pharmacogenomic approach." Bioinformatics 28.18 (2012): i611-i618.
	''')
	# get the drugs classified as PrimarySuspect or SecondarySuspect
	all_drugs = list(map(lambda entry: parse_drugs(entry), tqdm(all_drugs)))
	all_drugs = [x for x in tqdm(all_drugs) if x!=None]
	logging.info(f'Found {len(set([data[1] for data in all_drugs]))} unique drugs with PS or SS role')
	
	####################### BINDINGDB -- DB annotation specific ###########################
	# create the dicctionary with the drugs up to date for BINDINGDB-DrugBank
	annotation_file = f'./../../../Data/cross_side_information_DB/FDA/fda_bindingDB_db_dict{START_YEAR}_{END_YEAR}.txt'
	if os.path.isfile(annotation_file):
		logging.info(f'Found {annotation_file}')
		fda_DB_dict = read_bindingDB_annotation(annotation_file)
	else: 
		logging.info(f'Parsing the DB')
		tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
		root = tree.getroot()
		# db_drugProd_drugId
		bindingDB_dic = {} # dict with generic or product name as key and drugbank_id as value
		for drug_entry in tqdm(root):
			drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
			name = drug_entry.find('{http://www.drugbank.ca}name').text.upper()
			prod = drug_entry.find('{http://www.drugbank.ca}products')
			prod_names = set([brandname.find('{http://www.drugbank.ca}name').text.upper() for brandname in prod])
			if name:
				bindingDB_dic[name] = drugbank_ID
			if len(prod_names) >= 1:
				for prod in prod_names:
					bindingDB_dic[prod] = drugbank_ID
		blacklist = ['TABLETS', 'CAPSULES', 'INJECTABLE', 'INJECTION', 'MG', 'VIAL', 'ML']
		fda_DB_dict = []
		for fda_id, drug_name in tqdm(all_drugs):
			if bindingDB_dic.get(drug_name, None):
				fda_DB_dict.append((fda_id, drug_name, bindingDB_dic.get(drug_name)))
			else:
				drug_name_cleaned = re.sub(r'\([^)]*\)', '', drug_name)
				drug_name_cleaned = re.sub("[(,)]","", drug_name_cleaned)
				drug_name_cleaned = re.sub("[\d]+","", drug_name_cleaned)
				drug_name_cleaned = [word for word in drug_name_cleaned.split() if word not in blacklist]
				drug_name_cleaned = ' '.join(drug_name_cleaned).strip()
				if bindingDB_dic.get(drug_name_cleaned, None):
					fda_DB_dict.append((fda_id, drug_name_cleaned, bindingDB_dic.get(drug_name_cleaned)))
		fda_DB_dict = list(set(fda_DB_dict))
		assert len(set([len(event) for event in fda_DB_dict])) == 1, 'The FDA annotation is not well formatted'
		with open(annotation_file, 'w') as f:
			_ =f.write('fda_incident, drugname, DB_id\n')
			for entry in fda_DB_dict:
				_ = f.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\n')

	############################################################################
	all_events =  get_events(START_YEAR, END_YEAR)

	# get the events for ONLY the kept drugs
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
	unique_keywords= set(freq_percentage.keys())
	DB_PATH = '/home/margaret/data/jfuente/DTI/Data/BindingDB/BindingDB_All_2021m11/BindingDB_Kd_filtered_columns_subset_HomoSapiens_noduplicities.tsv'
	drugs = get_drugs_bindingDB(DB_PATH)
	drugs = list(set(drugs))
	unique_keywords_vector = sorted(list(unique_keywords))

	#IdDrug dict
	fda_db_dict_drug_id = {}
	for idn, _, drugname in tqdm(fda_DB_dict):
		if drugname in fda_db_dict_drug_id:
			fda_db_dict_drug_id[drugname].append(idn)
		else:
			fda_db_dict_drug_id[drugname] = [idn]

	#IDEvents AERS dict
	id_events_aers = {}
	for idn, aersname in tqdm(all_events):
		if idn in id_events_aers:
			id_events_aers[idn].append(aersname)
		else:
			id_events_aers[idn] = [aersname]

	aers_freq = []
	for drug in tqdm(drugs):
		aers_freq.append(get_freqs(
								drug, 
								unique_keywords_vector, 
								fda_db_dict_drug_id, 
								freq_percentage, 
								id_events_aers))

	aers_bit = []
	for freq in aers_freq:
		aers_bit.append([1 if entry > 0 else 0 for entry in freq])

	path = check_and_create_folder(db_name)
	aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
	aers_freq_pd.to_csv(os.path.join(path, 'aers_freq.tsv'), sep='\t')
	aers_bit_pd =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector, index=drugs)
	aers_bit_pd.to_csv(os.path.join(path, 'aers_bit.tsv'), sep='\t')
	

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################


