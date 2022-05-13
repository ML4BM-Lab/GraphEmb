import os, sys, re
import logging
import zipfile
import urllib.request
from tqdm import tqdm
import requests
import itertools
import xml.etree.ElementTree as ET
import pubchempy as pcp
import multiprocessing as mp
from random import randint
from time import sleep
from parallelbar import progress_map


def get_DB_name(path):
    """
    This function returns the name of the DB.
    """
    DB_NAMES = [
        "bindingDB",
        "BindingDB",
        "Davis_et_al",
        "DrugBank",
        "E",
        "GPCR",
        "IC",
        "NR",
        "BIOSNAP",
    ]
    for db in DB_NAMES:
        if re.search(db, path):
            logging.info(f"Database: {db}")
            if db in ["E", "GPCR", "IC", "NR"]:
                db = os.path.join("Yamanashi_et_al_GoldStandard", db)
                return db
            else:
                return db
    logging.error(f"Database: {db} not found")
    sys.exit("Please provide a valid database")


def read_dtis(path):
    """ "
    Read the dti file and return the dti along with
    the numerical dictionary sorted by name
    """
    with open(path, "r") as f:
        dtis = f.readlines()
    dtis = [entry.strip().split("\t") for entry in dtis]
    all_elements = list(set(element for entry in dtis for element in entry))
    dictionary = {node: index for index, node in enumerate(sorted(all_elements))}
    return dtis, dictionary


def read_and_extract_drugs(path):
    drugs = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                drug, _ = line.split("\t")
                drugs.append(drug)
    return drugs


def get_drugs_bindingDB(path):
    with open(path, "r") as f:
        _ = next(f)
        drugs = f.readlines()
    drugs = [int(float(entry.strip().split("\t")[1])) for entry in drugs]
    return drugs


def get_drugs_Davis(path):
    logging.info(f"Parsing the DB")
    with open(path, "r") as fl:
        _ = next(fl)
        all_entries = fl.readlines()
    all_entries = [int(entry.split("\t")[1]) for entry in all_entries]
    return all_entries


def get_drugs_DrugBank(path):
    drugs = []
    with open(path, "r") as f:
        _ = next(f)
        drugs = f.readlines()
    drugs = [entry.strip().split("\t")[0] for entry in drugs]
    return drugs


def get_drugs_BIOSNAP(path):
    with open(path, "r") as f:
        _ = next(f)
        drugs = f.readlines()
    drugs = [entry.strip().split("\t")[0] for entry in drugs]
    return drugs


def dowload_FAERS(start_year=2004, end_year=2021):
    """
    Downloads the FAERS data from the FDA website
    """
    logging.info(f"Downloading FAERS data from {start_year} to {end_year}")
    # it was a change on 2012Q4 from aers to faers
    for year in range(start_year, end_year + 1):
        for quarter in range(1, 5):
            dets_file = f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}.zip"
            if not os.path.exists(dets_file):
                logging.info(f"Downloading FAERS data from {year}-Q{quarter}")
                url = f"https://fis.fda.gov/content/Exports/faers_ascii_{year}Q{quarter}.zip"
                logging.debug(f"Downloading {url}")
                try:
                    urllib.request.urlretrieve(url, dets_file)
                except urllib.error.HTTPError:
                    logging.warning(
                        f"falling back to aers instead of faers for {year}-Q{quarter}"
                    )
                    url = f"https://fis.fda.gov/content/Exports/aers_ascii_{year}Q{quarter}.zip"
                    urllib.request.urlretrieve(url, dets_file)
            else:
                logging.debug(f"Already existing file for {year}-Q{quarter}. Skipping")


def uncompress_FAERS(path_2_folder):
    """
    Uncompress the FAERS data
    """
    zipfiles = filter(lambda file: file.endswith("zip"), os.listdir(path_2_folder))
    for file in zipfiles:
        dest = os.path.join(path_2_folder, file.split(".")[0])
        if os.path.exists(dest):
            logging.debug(f"{dest} already exists. Skipping")
            continue
        os.mkdir(dest)
        logging.info(f"Uncompressing {file}")
        with zipfile.ZipFile(os.path.join(path_2_folder, file), "r") as zip_ref:
            zip_ref.extractall(dest)


def get_events(start_year=2004, end_year=2021):
    all_events = []
    for year in tqdm(range(start_year, end_year + 1)):
        for quarter in range(1, 5):
            logging.debug(f"Processing FAERS data from {year}-Q{quarter}")
            # someone screwed up the naming of the files on the FDA side....
            possible_files = [
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/REAC{str(year)[-2:]}Q{quarter}.txt",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/REAC{str(year)[-2:]}Q{quarter}.TXT",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/REAC{str(year)[-2:]}Q{quarter}.txt",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/REAC{str(year)[-2:]}Q{quarter}.TXT",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.TXT",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/reac{str(year)[-2:]}q{quarter}.txt",
            ]
            fl = [
                real_file for real_file in possible_files if os.path.isfile(real_file)
            ]
            if fl:
                fl = fl[0]
            else:
                logging.warning(f"No file found for {year}-Q{quarter}")
                continue
            logging.debug(f"Reading {fl}")
            with open(fl, "r") as f:
                _ = next(f)
                events = f.readlines()
            events_cleaned = []
            # later files from the fda, have 3 fields instead of 2
            for event in events:
                event = event.strip().split("$")
                if len(event) == 3:
                    events_cleaned.append([event[0], event[1].upper()])
                else:
                    events_cleaned.append([event[0], event[2].upper()])
            all_events.extend(events_cleaned)
            logging.debug(f"{len(all_events)} total events")
    return all_events


def get_drugs(start_year=2004, end_year=2021):
    all_drugs = []
    for year in tqdm(range(start_year, end_year + 1), desc="Reading FDA drugs"):
        for quarter in range(1, 5):
            logging.debug(f"Processing FAERS data from {year}-Q{quarter}")
            # someone screwed up the naming of the files on the FDA side....
            possible_files = [
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/DRUG{str(year)[-2:]}Q{quarter}.txt",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ASCII/DRUG{str(year)[-2:]}Q{quarter}.TXT",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/DRUG{str(year)[-2:]}Q{quarter}.txt",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/DRUG{str(year)[-2:]}Q{quarter}.TXT",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/drug{str(year)[-2:]}q{quarter}.TXT",
                f"./../../DB/Data/cross_side_information_DB/FDA/faers_ascii_{year}Q{quarter}/ascii/drug{str(year)[-2:]}q{quarter}.txt",
            ]
            fl = [
                real_file for real_file in possible_files if os.path.isfile(real_file)
            ]
            if fl:
                fl = fl[0]
            else:
                logging.warning(f"No file found for {year}-Q{quarter}")
                continue
            with open(fl, "r", encoding="ISO-8859-1") as f:
                _ = next(f)
                drugs = f.readlines()
            all_drugs.extend(drugs)
    return all_drugs


def read_annotation(annotation_file):
    with open(annotation_file, "r") as f:
        _ = next(f)
        dictionary = f.readlines()
    dictionary = [entry.strip().split("\t") for entry in dictionary]
    return dictionary


def parse_drugs(entry):
    drugs_to_keep = ["PS", "SS"]
    entry = entry.strip().split("$")
    try:
        if entry[2] in drugs_to_keep:
            return (entry[0], entry[3])
        elif entry[3] in drugs_to_keep:
            return (entry[0], entry[4])
    except IndexError:
        logging.warning(f"{entry} is not a valid drug entry")


def get_drugs_FDA(START_YEAR, END_YEAR):
    all_drugs = get_drugs(START_YEAR, END_YEAR)
    # unique_drugs = len(set([entry.split('$')[3] for entry in all_drugs]))
    # logging.debug(f'Found {unique_drugs} unique drugs Vs the 291.997 claimed on the paper. ACCEPTABLE')
    logging.info(
        """
		In this study, we used the drugs (compound names or product names) 
		labeled as PS or SS and searched the equivalent drugs registered in the KEGG database.
		CITE:
		Takarabe, Masataka, et al. "Drug target prediction using adverse event report systems: 
		a pharmacogenomic approach." Bioinformatics 28.18 (2012): i611-i618.
	"""
    )
    # get the drugs classified as PrimarySuspect or SecondarySuspect
    all_drugs = [parse_drugs(drug) for drug in tqdm(all_drugs, desc="Parsing drugs")]
    all_drugs = [x for x in tqdm(all_drugs, desc="Removing empty entries") if x != None]
    logging.info(
        f"Found {len(set([data[1] for data in all_drugs]))} unique drugs with PS or SS role"
    )
    return all_drugs

def check_and_create_folder(db_name,model_name):
    folder = os.path.join("./"+model_name+"/Data", db_name)
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder


def get_freqs(drug, keywords, fda_dict, keyword_freq_percent, id_aers_dict):
    drug_kw_vector = [0] * len(keywords)
    all_cases = fda_dict.get(drug)
    if not all_cases:
        return drug_kw_vector
    all_drug_events = [id_aers_dict.get(case, [0]) for case in set(all_cases)]
    all_drug_events = [x for x in all_drug_events if x != 0]
    all_drug_events = [item for sublist in all_drug_events for item in sublist]
    keywords_set = {x: i for i, x in enumerate(keywords)}
    for event in all_drug_events:
        if event in keywords_set:
            drug_kw_vector[keywords_set.get(event)] = keyword_freq_percent.get(event)
    return drug_kw_vector

def clean_name(blacklist, drug_name):
    drug_name_cleaned = re.sub(r"\([^)]*\)", "", drug_name)
    drug_name_cleaned = re.sub("[(,)]", "", drug_name_cleaned)
    drug_name_cleaned = re.sub("[\d]+", "", drug_name_cleaned)
    drug_name_cleaned = [
        word for word in drug_name_cleaned.split() if word not in blacklist
    ]
    drug_name_cleaned = " ".join(drug_name_cleaned).strip()
    return drug_name_cleaned


def get_frequencies(
    fda_DB_dict, all_events, freq_percentage, drugs, unique_keywords_vector
):
    # IdDrug dict
    fda_db_dict_drug_id = {}
    for idn, _, drugname in tqdm(fda_DB_dict, desc="Creating drug id dict"):
        if drugname in fda_db_dict_drug_id:
            fda_db_dict_drug_id[drugname].append(idn)
        else:
            fda_db_dict_drug_id[drugname] = [idn]

    # IDEvents AERS dict
    id_events_aers = {}
    for idn, aersname in tqdm(all_events, desc="Creating id events aers dict"):
        if idn in id_events_aers:
            id_events_aers[idn].append(aersname)
        else:
            id_events_aers[idn] = [aersname]

    aers_freq = []
    for drug in tqdm(drugs, desc="Getting frequencies"):
        aers_freq.append(
            get_freqs(
                drug,
                unique_keywords_vector,
                fda_db_dict_drug_id,
                freq_percentage,
                id_events_aers,
            )
        )

    aers_bit = []
    for freq in aers_freq:
        aers_bit.append([1 if entry > 0 else 0 for entry in freq])
    return aers_freq, aers_bit


def read_and_extract_drugs_Yamanashi(path):
    drugs = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                drug, _ = line.split("\t")
                drugs.append(drug)
    return drugs


def get_dict_SIDER(path):
    with open(path, "r") as f:
        dictionary = f.readlines()
    dictionary = [(entry.strip().split("\t")) for entry in dictionary]
    return dict(dictionary)


def get_Kegg_pubchem_dict():
    PUBCHEM_PATTERN = re.compile(r"(?<=pubchem:)[\d]+")
    KEGG_PATTERN = re.compile(r"(?<=dr:)[\w\d]+")
    r = requests.get(f"http://rest.kegg.jp/conv/drug/pubchem")
    if r.status_code == 200:
        dictionary = r.text.strip().split("\n")
        dictionary = [
            (KEGG_PATTERN.search(entry).group(), PUBCHEM_PATTERN.search(entry).group())
            for entry in dictionary
        ]
        return dict(dictionary)


def get_side_effects_sider(path):
    CID_PATTERN = re.compile(r"(?<=CID1)[\d]+")
    with open(path, "r") as f:
        dictionary = f.readlines()
    dictionary = [(entry.strip().split("\t")) for entry in dictionary]
    dictionary = [(entry[0], entry[-1]) for entry in dictionary if entry[3] == "PT"]
    sider_SE_event_dict = {}
    for id, se in tqdm(dictionary, desc="Creating dict ID-SE"):
        id = CID_PATTERN.search(id).group()
        if id in sider_SE_event_dict:
            sider_SE_event_dict[id].append(se)
        else:
            sider_SE_event_dict[id] = [se]
    return sider_SE_event_dict


def get_Kegg_SIDER_dict(path):
    logging.info("Creating correlations KEGG-STITCH")
    CID_PATTERN = re.compile(r"(?<=CIDm)[\d]+")
    with open(path, "r") as fl:
        _ = next(fl)
        dictionary = fl.readlines()
    dictionary = [(entry.strip().split("\t")) for entry in dictionary]
    dictionary = [
        (entry[-1], CID_PATTERN.search(entry[0]).group()) for entry in dictionary
    ]
    return dict(dictionary)


def get_yamanashi_subDB(path):
    return re.search(r"(?<=\/)[\w]+", path).group()


def get_PubChem_SIDER_dict(path):
    logging.info("Creating correlations Davis-STITCH")
    CID_PATTERN = re.compile(r"(?<=CIDm)[\d]+")
    dictionary = {}
    with open(path, "r") as fl:
        _ = next(fl)
        for line in tqdm(fl, desc="Creating dictionary", total=70_000_000):
            line = line.split("\t")
            stitch = CID_PATTERN.search(line[0]).group()
            pc_id = line[-1].strip()
            dictionary[pc_id] = stitch
            # if pc_id in dictionary:
            # 	dictionary[pc_id].append(stitch)
            # else:
            # 	dictionary[pc_id] = [stitch]
    return dictionary


def get_Binding_to_Pubchem(drugs):
    tree = ET.parse(
        "./../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml"
    )
    root = tree.getroot()
    binding_to_pbC = {}
    for drug_entry in tqdm(root, desc="retrieving Pubchem IDs"):
        drugbank_ID = drug_entry.find("{http://www.drugbank.ca}drugbank-id").text
        if drugbank_ID in drugs:
            external_ids = drug_entry.find(
                "{http://www.drugbank.ca}external-identifiers"
            )
            if external_ids is not None:
                all_ext_ids = [ext_id[0].text for ext_id in external_ids]
            if "PubChem Compound" in all_ext_ids:
                binding_to_pbC[drugbank_ID] = external_ids[
                    all_ext_ids.index("PubChem Compound")
                ][1].text
                continue
            elif "PubChem Substance" in all_ext_ids:
                binding_to_pbC[drugbank_ID] = external_ids[
                    all_ext_ids.index("PubChem Substance")
                ][1].text
                continue
            else:
                logging.debug(
                    "Drug from DB not found in PubChem: {}".format(drugbank_ID)
                )
    return binding_to_pbC


def get_STITCH_from_Pubchem(compounds, substances):
    pubchem_to_stitch = {}
    with open(
        "./../../DB/Data/cross_side_information_DB/STITCH/dict_STITCH_PC.tsv", "r"
    ) as fl:
        _ = next(fl)
        for line in tqdm(
            fl, desc="Creating dictionary with Compounds", total=68_373_242
        ):
            ln = line.strip().split("\t")
            if ln[-1] in compounds:
                pubchem_to_stitch[ln[-1]] = ln[0]

    with open(
        "./../../DB/Data/cross_side_information_DB/STITCH/dict_STITCH_PS.tsv", "r"
    ) as fl:
        _ = next(fl)
        for line in tqdm(
            fl, desc="Creating dictionary with Substances", total=192_713_372
        ):
            ln = line.strip().split("\t")
            if ln[-1] in substances:
                pubchem_to_stitch[ln[-1]] = ln[0]

    return pubchem_to_stitch


# def get_names_from_Pubchem(drugs):
# 	synonyms = 3
# 	pbchem_to_names = {}
# 	for drug in tqdm(drugs, desc='Retrieving names from Pubchem'):
# 		names = pcp.get_synonyms(str(drug))
# 		sleep(randint(1,5))
# 		if names:
# 			names = names[0].get('Synonym')[:synonyms]
# 			for name in names:
# 				pbchem_to_names[name] = drug
# 		else:
# 			id = pcp.Substance.from_sid(drug).cids
# 			if  id:
# 				names = pcp.get_synonyms(id)
# 				if names:
# 					names = names[0].get('Synonym')[:synonyms]
# 					for name in names:
# 						pbchem_to_names[name] = drug
# 				else:
# 					logging.debug('Drug not found in Pubchem: {}'.format(drug))
# 			else:
# 					logging.debug('Drug not found in Pubchem: {}'.format(drug))
# 	return pbchem_to_names


# def get_names_from_Pubchem_mp(drugs, ncores):
# 	results = progress_map(get_names_from_pubchem, drugs, n_cpu=ncores, chunk_size=1, core_progress=True)
# 	# with mp.Pool(ncores) as pool:
# 	# 	results = pool.map(get_names_from_pubchem, drugs)
# 	return results


# def get_names_from_pubchem(drug):
# 	synonyms = 3
# 	names = pcp.get_synonyms(str(drug))
# 	sleep(randint(1,5))
# 	if names:
# 		names = names[0].get('Synonym')[:synonyms]
# 		return((drug, names))
# 	else:
# 		id = pcp.Substance.from_sid(drug).cids
# 		if  id:
# 			names = pcp.get_synonyms(id)
# 			if names:
# 				names = names[0].get('Synonym')[:synonyms]
# 				return((drug, names))


def getDB_from_binding(drugs, root=None):
    binding_to_DB = {}
    if not root:
        tree = ET.parse(
            "./../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml"
        )
        root = tree.getroot()
    for drug_entry in tqdm(root, desc="Retrieving DB IDs"):
        drugbank_ID = drug_entry.find("{http://www.drugbank.ca}drugbank-id").text
        external_ids = drug_entry.find("{http://www.drugbank.ca}external-identifiers")
        if external_ids is not None:
            all_ext_ids = [ext_id[0].text for ext_id in external_ids]
        if "PubChem Compound" in all_ext_ids:
            pubchem_compound = int(
                external_ids[all_ext_ids.index("PubChem Compound")][1].text
            )
            if pubchem_compound in drugs:
                binding_to_DB[pubchem_compound] = drugbank_ID
                continue
        elif "PubChem Substance" in all_ext_ids:
            pubchem_substance = external_ids[all_ext_ids.index("PubChem Substance")][
                1
            ].text
            if pubchem_substance in drugs:
                binding_to_DB[pubchem_substance] = drugbank_ID
                continue
    return binding_to_DB


def get_drugNames_from_Pubchem_web(drugs):
    puchchem_id_name = {}
    SYNONYMS = 5
    # split the list in chunks of 1000 to make requests
    size = 100
    drugs = [str(drug) for drug in drugs]
    split_list = lambda big_list, x: [
        big_list[i : i + x] for i in range(0, len(big_list), x)
    ]
    drug_chunks = split_list(drugs, size)
    for chunk in tqdm(drug_chunks):
        chunk = ",".join(chunk)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{chunk}/synonyms/json"
        response = requests.get(url)
        if response.status_code == 200:
            jsons = response.json()
        if response.status_code == 404:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/sid/{chunk}/synonyms/json"
            response = requests.get(url)
            jsons = response.json()
        for id in jsons.get("InformationList").get("Information"):
            cid = id.get("CID")
            names = id.get("Synonym")
            if names:
                names = names[:SYNONYMS]
                if cid in puchchem_id_name:
                    puchchem_id_name[cid].append(names)
                else:
                    puchchem_id_name[cid] = names
    return puchchem_id_name


def get_drugNames_from_Pubchem_web_inv(drugs):
    puchchem_id_name = {}
    # SYNONYMS = 5
    # split the list in chunks of 1000 to make requests
    size = 100
    drugs = [str(drug) for drug in drugs]
    split_list = lambda big_list, x: [
        big_list[i : i + x] for i in range(0, len(big_list), x)
    ]
    drug_chunks = split_list(drugs, size)
    for chunk in tqdm(drug_chunks):
        chunk = ",".join(chunk)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{chunk}/synonyms/json"
        response = requests.get(url)
        if response.status_code == 200:
            jsons = response.json()
        elif response.status_code == 404:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/sid/{chunk}/synonyms/json"
            response = requests.get(url)
            jsons = response.json()
        for id in jsons.get("InformationList").get("Information"):
            cid = id.get("CID")
            names = id.get("Synonym")
            if names:
                # names = names[:SYNONYMS]
                for name in names:
                    if name in puchchem_id_name:
                        puchchem_id_name[name].append(cid)
                    else:
                        puchchem_id_name[name] = [cid]
    return puchchem_id_name