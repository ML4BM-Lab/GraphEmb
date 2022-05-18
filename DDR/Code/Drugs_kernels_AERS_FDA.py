import os, sys
import re
import logging
import requests as rq
import argparse
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
from collections import Counter
import xml.etree.ElementTree as ET
import Drugs_kernels_helpers as hf


def DAVIS_AERS_FDA(DB_PATH = './DB/Data/Davis_et_al/Davis_30.tsv',model_name='DDR'):
	
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	# sanity check for the DB
	paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
	logging.info(f'\n{paper_cite}\n')

	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	START_YEAR = 2004
	END_YEAR   = 2021
	hf.dowload_FAERS(START_YEAR,END_YEAR)
	hf.uncompress_FAERS('./DB/Data/cross_side_information_DB/FDA')
	all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
	
	####################### DAVIS -- DB annotation specific ###########################
	# create the dicctionary with the drugs up to date for BINDINGDB-DrugBank
	annotation_file = f'./DB/Data/cross_side_information_DB/FDA/Davis_Drug_FDA_dict{START_YEAR}_{END_YEAR}.txt'
	if os.path.isfile(annotation_file):
		logging.info(f'Found {annotation_file}')
		fda_DB_dict = hf.read_annotation(annotation_file)
	else: 
		all_entries = hf.get_drugs_Davis(DB_PATH)
		# get the comercial names from pubchem
		NUMBER_OF_SYNONYMS = 3
		davis_PC_dic = {}
		for entry in tqdm(set(all_entries)):
			names = pcp.get_synonyms(entry)
			if names:
				names = names[0].get('Synonym')
				for name in names[:NUMBER_OF_SYNONYMS]:
					davis_PC_dic[name] = entry
		blacklist = ['TABLETS', 'CAPSULES', 'INJECTABLE', 'INJECTION', 'MG', 'VIAL', 'ML']
		fda_DB_dict = []
		for fda_id, drug_name in tqdm(all_drugs, desc='Crossing IDs with the FDA'):
			if davis_PC_dic.get(drug_name, None):
				fda_DB_dict.append((fda_id, drug_name, davis_PC_dic.get(drug_name)))
			else:
				drug_name_cleaned = hf.clean_name(blacklist, drug_name)
				if davis_PC_dic.get(drug_name_cleaned, None):
					fda_DB_dict.append((fda_id, drug_name_cleaned, davis_PC_dic.get(drug_name_cleaned)))
		
		fda_DB_dict = list(set(fda_DB_dict))
		assert len(set([len(event) for event in fda_DB_dict])) == 1, 'The FDA annotation is not well formatted'
		with open(annotation_file, 'w') as f:
			_ =f.write('fda_incident, drugname, DB_id\n')
			for entry in fda_DB_dict:
				_ = f.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\n')

	############################################################################
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
	
	drugs =  hf.get_drugs_Davis(DB_PATH)
	drugs = list(set(drugs))
	unique_keywords_vector = sorted(list(unique_keywords))

	aers_freq, aers_bit = hf.get_frequencies(fda_DB_dict, 
											all_events, 
											freq_percentage, 
											drugs, 
											unique_keywords_vector)

	path = hf.check_and_create_folder(db_name,model_name)
	aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
	aers_freq_pd.to_csv(os.path.join(path, 'Davis_et_al_drug_FDA_aers_freq.tsv'), sep='\t')
	aers_bit_pd =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector, index=drugs)
	aers_bit_pd.to_csv(os.path.join(path, 'Davis_et_al_drug_FDA_aers_bit.tsv'), sep='\t')

def BINDINGDB_AERS_FDA(DB_PATH = './DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv',model_name='DDR'):

    print(f'Actual path is {os.getcwd()}')

    fmt = '[%(levelname)s] %(message)s'
    logging.basicConfig(format=fmt, level=logging.DEBUG)

    # sanity check for the DB
    paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
    logging.info(f'\n{paper_cite}\n')
 
    logging.info(f'Reading database from: {DB_PATH}')
    db_name = hf.get_DB_name(DB_PATH)
    
    START_YEAR = 2004
    END_YEAR = 2021

    hf.dowload_FAERS(START_YEAR,END_YEAR)
    hf.uncompress_FAERS('./DB/Data/cross_side_information_DB/FDA')

    all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
    
    ####################### BINDINGDB -- DB annotation specific ###########################
    # create the dicctionary with the drugs up to date for BINDINGDB-DrugBank
    annotation_file = f'./DB/Data/cross_side_information_DB/FDA/fda_bindingDB_db_dict{START_YEAR}_{END_YEAR}.txt'
    if os.path.isfile(annotation_file):
        logging.info(f'Found {annotation_file}')
        fda_DB_dict = hf.read_annotation(annotation_file)
    else:
        db_puchem_dict_inv = hf.get_drugNames_from_Pubchem_web_inv(drugs)
        blacklist = ['TABLETS', 'CAPSULES', 'INJECTABLE', 'INJECTION', 'MG', 'VIAL', 'ML']
        fda_DB_dict = []
        for fda_id, drug_name in tqdm(all_drugs, desc = 'Getting names from FDA to Pubchem'):
            if db_puchem_dict_inv.get(drug_name.upper(), None):
                fda_DB_dict.append((fda_id, drug_name, db_puchem_dict_inv.get(drug_name.upper())[0]))
            else:
                drug_name_cleaned = re.sub(r'\([^)]*\)', '', drug_name)
                drug_name_cleaned = re.sub("[(,)]","", drug_name_cleaned)
                drug_name_cleaned = re.sub("[\d]+","", drug_name_cleaned)
                drug_name_cleaned = [word for word in drug_name_cleaned.split() if word not in blacklist]
                drug_name_cleaned = ' '.join(drug_name_cleaned).strip()
                if db_puchem_dict_inv.get(drug_name_cleaned.upper(), None):
                    fda_DB_dict.append((fda_id, drug_name_cleaned, db_puchem_dict_inv.get(drug_name_cleaned.upper())[0]))

        fda_DB_dict = list(set(fda_DB_dict))
        assert len(set([len(event) for event in fda_DB_dict])) == 1, 'The FDA annotation is not well formatted'
        with open(annotation_file, 'w') as f:
            _ =f.write('fda_incident, drugname, DB_id\n')
            for entry in fda_DB_dict:
                _ = f.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\n')

    ############################################################################
    all_events =  hf.get_events(START_YEAR, END_YEAR)

    # get the events for ONLY the kept drugs
    ids_2_keep = set([drug_entry[0]  for drug_entry in all_drugs])
    all_events = [event for event in tqdm(all_events, desc= 'Keeping the events for the selected drugs') if event[0] in ids_2_keep]

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
    drugs = hf.get_drugs_bindingDB(DB_PATH)
    drugs = list(set(drugs))
    # get the db_id for the pubchem_IDs drugs  # pickling error
    # pubchem_2_db = hf.get_names_from_Pubchem_mp(drugs, 20) # pickling error
    # pubchem_2_db = hf.get_names_from_Pubchem(drugs) # urllib.error.URLError: <urlopen error EOF occurred in violation of protocol (_ssl.c:1131)>
    # db_puchem_dict = hf.get_drugNames_from_Pubchem_web(drugs)
    db_puchem_dict_inv = hf.get_drugNames_from_Pubchem_web_inv(drugs)
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

            # binding_2_pc = hf.get_Binding_to_Pubchem(drugs)
    aers_freq = []
    for drug in tqdm(drugs):
        aers_freq.append(hf.get_freqs(
                                drug, 
                                unique_keywords_vector, 
                                fda_db_dict_drug_id, 
                                freq_percentage, 
                                id_events_aers))
    
    if sum([sum(entry) for entry in aers_freq]) != 0:
        logging.warning('The AERS frequencies are all zeroes.')

    aers_bit = []
    for freq in aers_freq:
        aers_bit.append([1 if entry > 0 else 0 for entry in freq])

    path = hf.check_and_create_folder(db_name,model_name)
    logging.info(f'Saving the results in {path}')
    aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
    aers_freq_pd.to_csv(os.path.join(path, 'BindingDB_drug_FDA_aers_freq.tsv'), sep='\t')
    #aers_freq_pd.to_pickle(os.path.join(path, 'BindingDB_Drug_FDA_aers_freq.pickle'))
    aers_bit_pd =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector, index=drugs)
    aers_bit_pd.to_csv(os.path.join(path, 'BindingDB_drug_FDA_aers_bit.tsv'), sep='\t')
    #aers_bit_pd.to_pickle(os.path.join(path, 'BindingDB_Drug_FDA_aers_bit.pickle'))

def DRUGBANK_AERS_FDA(DB_PATH =  './DB/Data/DrugBank/DrugBank_DTIs.tsv',model_name='DDR'):
	
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	# sanity check for the DB
	paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
	logging.info(f'\n{paper_cite}\n')
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	START_YEAR = 2004
	END_YEAR   = 2021
	hf.dowload_FAERS(START_YEAR,END_YEAR)
	hf.uncompress_FAERS('./DB/Data/cross_side_information_DB/FDA')
	all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
	
	####################### DRUGBANK -- DB annotation specific ###########################
	# create the dicctionary with the drugs up to date for BINDINGDB-DrugBank
	annotation_file = f'./DB/Data/cross_side_information_DB/FDA/DrugBank_Drug_FDA_dict{START_YEAR}_{END_YEAR}.txt'
	if os.path.isfile(annotation_file):
		logging.info(f'Found {annotation_file}')
		fda_DB_dict = hf.read_annotation(annotation_file)
	else: 
		all_entries = list(set(hf.get_drugs_DrugBank(DB_PATH)))
		# get the names from drugbank
		tree = ET.parse('./DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
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

	path = hf.check_and_create_folder(db_name,model_name)
	aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
	aers_freq_pd.to_csv(os.path.join(path, 'DrugBank_drug_FDA_aers_freq.tsv'), sep='\t')
	aers_bit_pd =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector, index=drugs)
	aers_bit_pd.to_csv(os.path.join(path, 'DrugBank_drug_FDA_aers_bit.tsv'), sep='\t')

def YAMANASHI_AERS_FDA(subdataset='E',model_name='DDR'):
    
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)
	DB_PATH = './DB/Data/Yamanashi_et_al_GoldStandard/'+subdataset+'/interactions/'+subdataset.lower()+'_admat_dgc_mat_2_line.txt'

	# sanity check for the DB
	paper_cite = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/'
	logging.info(f'\n\n{paper_cite}\n')
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)
		
	START_YEAR = 2004
	END_YEAR = 2021

	hf.dowload_FAERS(START_YEAR,END_YEAR)
	hf.uncompress_FAERS('./DB/Data/cross_side_information_DB/FDA')

	all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
	
	####################### Yamanashi -- KEGG annotation specific ###########################
	# create the dicctionary with the drugs up to date for KEGG
	annotation_file = f'./DB/Data/cross_side_information_DB/FDA//fda_kegg_dict{START_YEAR}_{END_YEAR}.txt'
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
	path = hf.check_and_create_folder(db_name,model_name)
	subname = hf.get_yamanashi_subDB(db_name)
	aers_freq.to_csv(os.path.join(path, f'{subname}_drug_FDA_aers_freq.tsv'), sep='\t')
	aers_bit.to_csv(os.path.join(path,  f'{subname}_drug_FDA_aers_bit.tsv'), sep='\t')

def BIOSNAP_AERS_FDA(DB_PATH = '/home/margaret/data/jfuente/DTI/Input4Models/DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv',model_name='DDR'):

	
    fmt = "[%(levelname)s] %(message)s"
    logging.basicConfig(format=fmt, level=logging.DEBUG)

    # sanity check for the DB
    paper_cite = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436840/"
    logging.info(f"\n{paper_cite}\n")
    logging.info(f"Reading database from: {DB_PATH}")
    db_name = hf.get_DB_name(DB_PATH)

    START_YEAR = 2004
    END_YEAR = 2021

    hf.dowload_FAERS(START_YEAR, END_YEAR)
    hf.uncompress_FAERS("./DB/Data/cross_side_information_DB/FDA")
    all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)

    ####################### BIOSNAP -- DB annotation specific ###########################
    # create the dicctionary with the drugs up to date for BIOSNAP-DrugBank
    annotation_file = f"./DB/Data/cross_side_information_DB/FDA/fda_biosnap_db_dict{START_YEAR}_{END_YEAR}.txt"
    if os.path.isfile(annotation_file):
        logging.info(f"Found {annotation_file}")
        fda_DB_dict = hf.read_annotation(annotation_file)
    else:
        logging.info(f"Parsing the DB")
        tree = ET.parse(
            "./DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml"
        )
        root = tree.getroot()
        # db_drugProd_drugId
        biosnap_dic = (
            {}
        )  # dict with generic or product name as key and drugbank_id as value
        for drug_entry in tqdm(root, desc="Parsing the Drugbank"):
            drugbank_ID = drug_entry.find("{http://www.drugbank.ca}drugbank-id").text
            name = drug_entry.find("{http://www.drugbank.ca}name").text.upper()
            prod = drug_entry.find("{http://www.drugbank.ca}products")
            prod_names = set(
                [
                    brandname.find("{http://www.drugbank.ca}name").text.upper()
                    for brandname in prod
                ]
            )
            if name:
                biosnap_dic[name] = drugbank_ID
            if len(prod_names) >= 1:
                for prod in prod_names:
                    biosnap_dic[prod] = drugbank_ID
        blacklist = [
            "TABLETS",
            "CAPSULES",
            "INJECTABLE",
            "INJECTION",
            "MG",
            "VIAL",
            "ML",
        ]
        fda_DB_dict = []
        for fda_id, drug_name in tqdm(all_drugs, desc="Merging the databases"):
            if biosnap_dic.get(drug_name, None):
                fda_DB_dict.append((fda_id, drug_name, biosnap_dic.get(drug_name)))
            else:
                drug_name_cleaned = re.sub(r"\([^)]*\)", "", drug_name)
                drug_name_cleaned = re.sub("[(,)]", "", drug_name_cleaned)
                drug_name_cleaned = re.sub("[\d]+", "", drug_name_cleaned)
                drug_name_cleaned = [
                    word for word in drug_name_cleaned.split() if word not in blacklist
                ]
                drug_name_cleaned = " ".join(drug_name_cleaned).strip()
                if biosnap_dic.get(drug_name_cleaned, None):
                    fda_DB_dict.append(
                        (fda_id, drug_name_cleaned, biosnap_dic.get(drug_name_cleaned))
                    )

        fda_DB_dict = list(set(fda_DB_dict))
        assert (
            len(set([len(event) for event in fda_DB_dict])) == 1
        ), "The FDA annotation is not well formatted"
        with open(annotation_file, "w") as f:
            _ = f.write("fda_incident, drugname, DB_id\n")
            for entry in tqdm(fda_DB_dict, desc="Writing the annotation file"):
                _ = f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\n")

    ############################################################################
    all_events = hf.get_events(START_YEAR, END_YEAR)

    # get the events for ONLY the kept drugs
    ids_2_keep = set([drug_entry[0] for drug_entry in all_drugs])
    all_events = [event for event in tqdm(all_events) if event[0] in ids_2_keep]

    # get the keywords from the events
    all_keywords = [element for event in all_events for element in event]
    all_keywords = [
        element for element in tqdm(all_keywords) if bool(re.search(r"\D", element))
    ]

    freq = Counter(all_keywords)
    logging.info(
        """ Keywords that appear too frequently (>0.001) or too rarely (the number of reports were <5) were removed.
		remove the reports were <5
		CITE: 
		Takarabe, Masataka, et al. "Drug target prediction using adverse event report systems: 
		a pharmacogenomic approach." Bioinformatics 28.18 (2012): i611-i618.
		"""
    )
    freq_percentage = {
        drug_id: envents / len(all_keywords)
        for drug_id, envents in freq.items()
        if envents > 5
    }
    freq_percentage = {
        drug_id: envents
        for drug_id, envents in freq_percentage.items()
        if envents < 0.001
    }
    unique_keywords = set(freq_percentage.keys())

    drugs = hf.get_drugs_BIOSNAP(DB_PATH)
    drugs = list(set(drugs))
    unique_keywords_vector = sorted(list(unique_keywords))

    aers_freq, aers_bit = hf.get_frequencies(
        fda_DB_dict, all_events, freq_percentage, drugs, unique_keywords_vector
    )

    path = hf.check_and_create_folder(db_name,model_name)
    aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
    aers_freq_pd.to_csv(os.path.join(path, "BIOSNAP_drug_FDA_aers_freq.tsv"), sep="\t")
    aers_bit_pd = pd.DataFrame(aers_bit, columns=unique_keywords_vector, index=drugs)
    aers_bit_pd.to_csv(os.path.join(path, "BIOSNAP__drug_FDA_aers_bit.tsv"), sep="\t")

