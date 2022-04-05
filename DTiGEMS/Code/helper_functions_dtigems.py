
import os, sys, re
import logging
import zipfile
import urllib.request
from tqdm import tqdm
import requests

def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['bindingDB', 'BindingDB', 'Davis_et_al', 'DrugBank', 'E', 'GPCR', 'IC', 'NR']
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

def get_drugs_Davis(path):
	logging.info(f'Parsing the DB')
	with open(path, 'r') as fl:
		_ = next(fl)
		all_entries = fl.readlines()
	all_entries = [int(entry.split('\t')[1]) for entry in all_entries]
	return all_entries

def get_drugs_DrugBank(path):
	drugs = []
	with open(path, 'r') as f:
		_ = next(f)
		drugs = f.readlines()
	drugs = [entry.strip().split('\t')[0] for entry in drugs]
	return drugs

def dowload_FAERS( start_year = 2004, end_year = 2021):
	"""
	Downloads the FAERS data from the FDA website
	"""
	logging.info(f'Downloading FAERS data from {start_year} to {end_year}')
	# it was a change on 2012Q4 from aers to faers
	for year in range( start_year, end_year+1):
		for quarter in range(1, 5):
			dets_file = f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}.zip'
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
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/REAC{str(year)[-2:]}Q{quarter}.txt',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/REAC{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/REAC{str(year)[-2:]}Q{quarter}.txt',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/REAC{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.TXT',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.txt'
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
	for year in tqdm(range(start_year,end_year+1), desc='Reading FDA drugs'):
		for quarter in range(1, 5):
			logging.debug(f'Processing FAERS data from {year}-Q{quarter}')
			# someone screwed up the naming of the files on the FDA side....
			possible_files = [
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/DRUG{str(year)[-2:]}Q{quarter}.txt',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/DRUG{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/DRUG{str(year)[-2:]}Q{quarter}.txt',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/DRUG{str(year)[-2:]}Q{quarter}.TXT',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/drug{str(year)[-2:]}q{quarter}.TXT',
				f'./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/drug{str(year)[-2:]}q{quarter}.txt'
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

def read_annotation(annotation_file):
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

def get_drugs_FDA(START_YEAR, END_YEAR):
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
	all_drugs = [parse_drugs(drug) for drug in tqdm(all_drugs, desc='Parsing drugs')]
	all_drugs = [x for x in tqdm(all_drugs, desc='Removing empty entries') if x!=None]
	logging.info(f'Found {len(set([data[1] for data in all_drugs]))} unique drugs with PS or SS role')
	return all_drugs

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

def clean_name(blacklist, drug_name):
    drug_name_cleaned = re.sub(r'\([^)]*\)', '', drug_name)
    drug_name_cleaned = re.sub("[(,)]","", drug_name_cleaned)
    drug_name_cleaned = re.sub("[\d]+","", drug_name_cleaned)
    drug_name_cleaned = [word for word in drug_name_cleaned.split() if word not in blacklist]
    drug_name_cleaned = ' '.join(drug_name_cleaned).strip()
    return drug_name_cleaned

def get_frequencies(fda_DB_dict, all_events, freq_percentage, drugs, unique_keywords_vector):
	#IdDrug dict
	fda_db_dict_drug_id = {}
	for idn, _, drugname in tqdm(fda_DB_dict, desc='Creating drug id dict'):
		if drugname in fda_db_dict_drug_id:
			fda_db_dict_drug_id[drugname].append(idn)
		else:
			fda_db_dict_drug_id[drugname] = [idn]
	
	#IDEvents AERS dict
	id_events_aers = {}
	for idn, aersname in tqdm(all_events, desc='Creating id events aers dict'):
		if idn in id_events_aers:
			id_events_aers[idn].append(aersname)
		else:
			id_events_aers[idn] = [aersname]

	aers_freq = []
	for drug in tqdm(drugs, desc='Getting frequencies'):
		aers_freq.append(get_freqs( drug, 
								unique_keywords_vector, 
								fda_db_dict_drug_id, 
								freq_percentage, 
								id_events_aers))

	aers_bit = []
	for freq in aers_freq:
		aers_bit.append([1 if entry > 0 else 0 for entry in freq])
	return aers_freq,aers_bit

def read_and_extract_drugs_Yamanashi(path):
	drugs = []
	with open(path, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				drug, _ = line.split('\t')
				drugs.append(drug)
	return drugs

def get_dict_SIDER(path):
	with open(path, 'r') as f:
		dictionary = f.readlines()
	dictionary = [(entry.strip().split('\t')) for entry in dictionary]
	return dict(dictionary)

def get_Kegg_pubchem_dict():
	PUBCHEM_PATTERN = re.compile(r'(?<=pubchem:)[\d]+')
	KEGG_PATTERN = re.compile(r'(?<=dr:)[\w\d]+')
	r = requests.get(f'http://rest.kegg.jp/conv/drug/pubchem')
	if r.status_code == 200:
		dictionary =  r.text.strip().split('\n')
		dictionary = [(KEGG_PATTERN.search(entry).group(), PUBCHEM_PATTERN.search(entry).group()) for entry in dictionary]
		return dict(dictionary)

def get_side_effects_sider(path):
	CID_PATTERN = re.compile(r'(?<=CID1)[\d]+')
	with open(path, 'r') as f:
		dictionary = f.readlines()
	dictionary = [(entry.strip().split('\t')) for entry in dictionary]
	dictionary = [(entry[0], entry[-1] )for entry in dictionary if entry[3] == 'PT']
	sider_SE_event_dict = {}
	for id, se in tqdm(dictionary, desc='Creating dict ID-SE'):
		id = CID_PATTERN.search(id).group()
		if id in sider_SE_event_dict:
			sider_SE_event_dict[id].append(se)
		else:
			sider_SE_event_dict[id] = [se]
	return sider_SE_event_dict

def get_Kegg_SIDER_dict(path):
	logging.info('Creating correlations KEGG-STITCH')
	CID_PATTERN = re.compile(r'(?<=CIDm)[\d]+')
	with open(path, 'r') as fl:
		_ = next(fl)
		dictionary = fl.readlines()
	dictionary = [(entry.strip().split('\t')) for entry in dictionary]
	dictionary = [(entry[-1], CID_PATTERN.search(entry[0]).group())for entry in dictionary ]
	return dict(dictionary)

def get_yamanashi_subDB(path):
	return re.search(r'(?<=\/)[\w]+', path).group()

# def get_Davis_SIDER_dict(path):
# 	logging.info('Creating correlations Davis-STITCH')
# 	CID_PATTERN = re.compile(r'(?<=CIDm)[\d]+')
# 	with open(path, 'r') as fl:
# 		_ = next(fl)
# 		dictionary = fl.readlines()
# 	dictionary = [(entry.strip().split('\t')) for entry in dictionary]
# 	dictionary = [(entry[-1], CID_PATTERN.search(entry[0]).group()) for entry in tqdm(dictionary, desc='Creating dictionary')  ]
# 	return dict(dictionary)

def get_Davis_SIDER_dict(path):
	logging.info('Creating correlations Davis-STITCH')
	CID_PATTERN = re.compile(r'(?<=CIDm)[\d]+')
	dictionary = {}
	with open(path, 'r') as fl:
		_ = next(fl)
		for line in tqdm(fl, desc='Creating dictionary', total=70_000_000):
			line = line.split('\t')
			stitch = CID_PATTERN.search(line[0]).group()
			pc_id = line[-1].strip()
			dictionary[pc_id] = stitch
			# if pc_id in dictionary:
			# 	dictionary[pc_id].append(stitch)
			# else:
			# 	dictionary[pc_id] = [stitch]
	return dictionary
