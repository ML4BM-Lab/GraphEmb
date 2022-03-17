
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


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR']
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
				logging.info(f'Already existing file for {year}-Q{quarter}. Skipping')

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
		with open(fl, 'r') as f:
			_ = next(f)
			events = f.readlines()
		events = [event.strip().split('$') for event in events]
		all_events.extend(events)
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
				f'./../../../Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.txt'
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

def read_kegg_annotation(annotation_file):
	with open(annotation_file, 'r') as f:
		_ = next(f)
		fda_kegg_dict = f.readlines()
	fda_kegg_dict = [entry.strip().split('\t') for entry in fda_kegg_dict]
	return fda_kegg_dict

def parse_drugs(entry):
	entry = entry.strip().split('$')
	if entry[2] == 'PS' or entry[2] == 'SS':
		return True
	else:
		return False

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
	all_drug_events = reduce(lambda x,y: x+y, map(lambda x: id_aers_dict.get(x,[0]), set(all_cases)))
	all_drug_events = [x for x in all_drug_events if x!=0]
	for event in all_drug_events:
		if event in keywords:
			drug_kw_vector[keywords.index(event)] = keyword_freq_percent.get(event)
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

	level= log_levels[args.verbosity]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)
	# parse the ascii ones not the xmls

	# set the logging info
	level= log_levels[args.verbosity]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)

	# sanity check for the DB
	paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
	logging.info(f'\n\n{paper_cite}\n')
	DB_PATH = args.dbPath
	# DB_PATH =  '/home/margaret/data/jfuente/DTI/Data/Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt'
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
	all_drugs = list(filter(lambda entry: parse_drugs(entry), tqdm(all_drugs)))
	all_drugs = [(data.split('$')[0], str(data.split('$')[3]).upper()) for data in tqdm(all_drugs)]
	logging.info(f'Found {len(set([data[1] for data in all_drugs]))} unique drugs with PS or SS role')
	
	####################### Yamanashi -- KEGG annotation specific ###########################
	# create the dicctionary with the drugs up to date for KEGG
	annotation_file = f'./../../../Data/cross_side_information_DB/FDA/fda_kegg_dict{START_YEAR}_{END_YEAR}.txt'
	if os.path.isfile(annotation_file):
		logging.info(f'Found {annotation_file}')
		fda_kegg_dict = read_kegg_annotation(annotation_file)
	else: 
		kegg_drug_names = rq.get('http://rest.kegg.jp/list/drug/').text.split('\n')
		kegg_drug_names = [entry.upper() for entry in kegg_drug_names]
		PATTERN_DB_KEGG = re.compile(r"(?<=DR:)[\w]+")
		fda_kegg_dict = []
		for fda_id, drugname in tqdm(all_drugs):
			for kegg_entry in kegg_drug_names:
				if drugname in kegg_entry:
					kegg_id = PATTERN_DB_KEGG.findall(kegg_entry)[0]
					fda_kegg_dict.append((fda_id, drugname, kegg_id))
					break
		with open(annotation_file, 'w') as f:
			_ =f.write('fda_incident, drugname, kegg_id\n')
			for entry in tqdm(fda_kegg_dict):
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

	drugs = read_and_extract_drugs(DB_PATH)
	drugs = list(set(drugs))

	unique_keywords_vector = sorted(list(unique_keywords))

	#IdDrug dict
	fda_kegg_dict_drug_id = {}
	for idn, _, drugname in tqdm(fda_kegg_dict):
		if drugname in fda_kegg_dict_drug_id:
			fda_kegg_dict_drug_id[drugname].append(idn)
		else:
			fda_kegg_dict_drug_id[drugname] = [idn]

	#IDEvents AERS dict
	id_events_aers = {}
	for idn, aersname, _ in tqdm(all_events):
		if idn in id_events_aers:
			id_events_aers[idn].append(aersname)
		else:
			id_events_aers[idn] = [aersname]

	aers_freq = []
	for drug in tqdm(drugs):
		aers_freq.append(get_freqs(drug, 
									unique_keywords_vector, 
									fda_kegg_dict_drug_id, 
									freq_percentage, id_events_aers))

	aers_bit = []
	for freq in aers_freq:
		aers_bit.append([1 if entry > 0 else 0 for entry in freq])

	aers_bit =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector,index=drugs)
	aers_freq = pd.DataFrame(aers_freq, columns=unique_keywords_vector,index=drugs)
	path = check_and_create_folder(db_name)
	aers_freq.to_csv(os.path.join(path, 'aers_freq.tsv'), sep='\t')
	aers_bit.to_csv(os.path.join(path, 'aers_bit.tsv'), sep='\t')


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################