
import os, sys
import re
import logging
import pandas as pd
import requests as rq
from tqdm import tqdm
from collections import Counter
import multiprocessing as mp
import argparse
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
    # DB_PATH =  './../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
    logging.info(f'Reading database from: {DB_PATH}')
    db_name = hf.get_DB_name(DB_PATH)
    
    START_YEAR = 2004
    END_YEAR = 2021

    hf.dowload_FAERS(START_YEAR,END_YEAR)
    hf.uncompress_FAERS('./../../DB/Data/cross_side_information_DB/FDA')

    all_drugs = hf.get_drugs_FDA(START_YEAR, END_YEAR)
    
    ####################### BINDINGDB -- DB annotation specific ###########################
    # create the dicctionary with the drugs up to date for BINDINGDB-DrugBank
    annotation_file = f'./../../DB/Data/cross_side_information_DB/FDA/fda_bindingDB_db_dict{START_YEAR}_{END_YEAR}.txt'
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

    path = hf.check_and_create_folder(db_name)
    logging.info(f'Saving the results in {path}')
    aers_freq_pd = pd.DataFrame(aers_freq, columns=unique_keywords_vector, index=drugs)
    aers_freq_pd.to_csv(os.path.join(path, 'BindingDB_Drug_FDA_aers_freq.tsv'), sep='\t')
    aers_freq_pd.to_pickle(os.path.join(path, 'BindingDB_Drug_FDA_aers_freq.pickle'))
    aers_bit_pd =  pd.DataFrame(aers_bit,  columns=unique_keywords_vector, index=drugs)
    aers_bit_pd.to_csv(os.path.join(path, 'BindingDB_Drug_FDA_aers_bit.tsv'), sep='\t')
    aers_bit_pd.to_pickle(os.path.join(path, 'BindingDB_Drug_FDA_aers_bit.pickle'))


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################




