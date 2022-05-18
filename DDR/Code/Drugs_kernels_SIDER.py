from sklearn.metrics.pairwise import cosine_similarity
import logging
import argparse
import pandas as pd
import numpy as np
import os, re
from tqdm import tqdm
import Drugs_kernels_helpers as hf


def BINDINGDB_SIDER(DB_PATH = './DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv', model_name = 'DDR'):

	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	# sanity check for the DB
	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3'
	logging.info(f'\n{paper_cite}\n')
	# 
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	drugs = list(set(hf.get_drugs_bindingDB(DB_PATH)))

	sider_SE_dict = hf.get_side_effects_sider('./DB/Data/cross_side_information_DB/SIDER/meddra_all_se.tsv')
	davis_SIDER_dict = hf.get_PubChem_SIDER_dict('./DB/Data/cross_side_information_DB/STITCH/dict_STITCH_PC.tsv')

	drug_side_effects = []
	all_side_effects = []
	for drug in drugs:
		binding_stitch_ID = davis_SIDER_dict.get(str(drug), None)
		if binding_stitch_ID:
			sider_SE = sider_SE_dict.get(binding_stitch_ID, None)
			if sider_SE:
				all_side_effects.append(sider_SE)
				drug_side_effects.append((drug, sider_SE))
			else:
				logging.debug(f'{drug} with no SE')	
		else:
			logging.debug(f'{drug} not found in the SIDER database')

	drug_side_effects = dict(drug_side_effects)

	all_side_effects_positions = [element for sublist in all_side_effects for element in sublist]
	all_side_effects_positions = sorted(list(set(all_side_effects_positions)))
	all_side_effects_positions = {v:k for k,v in enumerate(all_side_effects_positions)}

	sider_drugs_bin = []
	for drug in sorted(drugs):
		side_effects = drug_side_effects.get(drug)
		drug_se = [0]*len(all_side_effects_positions)
		if  not side_effects:
			sider_drugs_bin.append(drug_se)
		else:
			for se in side_effects:
				pos = all_side_effects_positions.get(se)
				drug_se[pos] = 1
			sider_drugs_bin.append(drug_se)

	sider_bit =  pd.DataFrame(sider_drugs_bin,  columns=list(all_side_effects_positions.keys()),index=sorted(drugs))

	sider_bit =  pd.DataFrame(cosine_similarity(sider_bit), columns=sorted(drugs), index=sorted(drugs))
	np.fill_diagonal(sider_bit.values, 1)

	path = hf.check_and_create_folder(db_name,model_name)
	# sider_bit.to_pickle(os.path.join(path, 'Binding_drug_SIDER_SideEffect.pickle'))
	sider_bit.to_csv(os.path.join(path, 'Binding_drug_SIDER_SideEffect.tsv'), sep='\t')
	
def DAVIS_SIDER(DB_PATH = './DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al.tsv', model_name = 'DDR'):

	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	# sanity check for the DB
	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3'
	logging.info(f'\n{paper_cite}\n')
	
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	drugs = list(set(hf.get_drugs_Davis(DB_PATH)))

	sider_SE_dict = hf.get_side_effects_sider('./DB/Data/cross_side_information_DB/SIDER/meddra_all_se.tsv')
	davis_SIDER_dict = hf.get_PubChem_SIDER_dict('./DB/Data/cross_side_information_DB/STITCH/dict_STITCH_PC.tsv')

	drug_side_effects = []
	all_side_effects = []
	for drug in drugs:
		davis_stitch_ID = davis_SIDER_dict.get(str(drug), None)
		if davis_stitch_ID:
			sider_SE = sider_SE_dict.get(davis_stitch_ID, None)
			if sider_SE:
				all_side_effects.append(sider_SE)
				drug_side_effects.append((drug, sider_SE))
			else:
				logging.debug(f'{drug} with no SE')	
		else:
			logging.debug(f'{drug} not found in the SIDER database')

	drug_side_effects = dict(drug_side_effects)

	all_side_effects_positions = [element for sublist in all_side_effects for element in sublist]
	all_side_effects_positions = sorted(list(set(all_side_effects_positions)))
	all_side_effects_positions = {v:k for k,v in enumerate(all_side_effects_positions)}

	sider_drugs_bin = []
	for drug in sorted(drugs):
		side_effects = drug_side_effects.get(drug)
		drug_se = [0]*len(all_side_effects_positions)
		if  not side_effects:
			sider_drugs_bin.append(drug_se)
		else:
			for se in side_effects:
				pos = all_side_effects_positions.get(se)
				drug_se[pos] = 1
			sider_drugs_bin.append(drug_se)

	sider_bit =  pd.DataFrame(sider_drugs_bin,  columns=list(all_side_effects_positions.keys()),index=sorted(drugs))

	sider_bit =  pd.DataFrame(cosine_similarity(sider_bit), columns=sorted(drugs), index=sorted(drugs))
	np.fill_diagonal(sider_bit.values, 1)

	path = hf.check_and_create_folder(db_name,model_name)
	# sider_bit.to_pickle(os.path.join(path, 'Davis_drug_SIDER_SideEffect.pickle'))
	sider_bit.to_csv(os.path.join(path, 'Davis_drug_SIDER_SideEffect.tsv'), sep='\t')
	
def BIOSNAP_SIDER(DB_PATH =  './DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv', model_name = 'DDR'):

	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	# sanity check for the DB
	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3'
	logging.info(f'\n{paper_cite}\n')

	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	drugs = list(set(hf.get_drugs_BIOSNAP(DB_PATH)))

	sider_SE_dict = hf.get_side_effects_sider('./DB/Data/cross_side_information_DB/SIDER/meddra_all_se.tsv')
	biosnap_to_pbC, biosnap_to_pbS = hf.get_Biosnap_to_Pubchem(drugs)
	pubchem_to_stitch = hf.get_STITCH_from_Pubchem(set(biosnap_to_pbC.values()), set(biosnap_to_pbS.values()))

	drug_side_effects = []
	all_side_effects = []
	CID_PATTERN = re.compile(r'(?<=CID[ms]{1})[\d]+')
	for drug in tqdm(drugs, desc='Assigning Side Effects to drugs'):
		pubchem_id = biosnap_to_pbC.get(drug,None) if drug in biosnap_to_pbC else biosnap_to_pbS.get(drug,None) 
		stitch = pubchem_to_stitch.get(pubchem_id, None)
		if stitch:
			stitch = CID_PATTERN.search(stitch).group()
			sider_SE = sider_SE_dict.get(stitch)
			if sider_SE:
				all_side_effects.append(sider_SE)
				drug_side_effects.append((drug, sider_SE))
			else:
				logging.debug(f'{drug} with no SE')	
		else:
			logging.debug(f'{drug} not found in the SIDER database')

	drug_side_effects = dict(drug_side_effects)

	all_side_effects_positions = [element for sublist in all_side_effects for element in sublist]
	all_side_effects_positions = sorted(list(set(all_side_effects_positions)))
	all_side_effects_positions = {v:k for k,v in enumerate(all_side_effects_positions)}

	sider_drugs_bin = []
	for drug in sorted(drugs):
		side_effects = drug_side_effects.get(drug)
		drug_se = [0]*len(all_side_effects_positions)
		if  not side_effects:
			sider_drugs_bin.append(drug_se)
		else:
			for se in side_effects:
				pos = all_side_effects_positions.get(se)
				drug_se[pos] = 1
			sider_drugs_bin.append(drug_se)

	sider_bit =  pd.DataFrame(sider_drugs_bin,  columns=list(all_side_effects_positions.keys()),index=sorted(drugs))
	sider_bit =  pd.DataFrame(cosine_similarity(sider_bit), columns=sorted(drugs), index=sorted(drugs))
	np.fill_diagonal(sider_bit.values, 1)

	path = hf.check_and_create_folder(db_name,model_name)
	# sider_bit.to_pickle(os.path.join(path, 'BIOSNAP_Drug_SIDER_SideEffect.pickle'))
	sider_bit.to_csv(os.path.join(path, 'BIOSNAP_drug_SIDER_SideEffect.tsv'), sep='\t')	

def DRUGBANK_SIDER(DB_PATH =  './DB/Data/DrugBank/DrugBank_DTIs.tsv', model_name = 'DDR'):

    # set the logging info
    fmt = "[%(levelname)s] %(message)s"
    logging.basicConfig(format=fmt, level=logging.DEBUG)

    # sanity check for the DB
    paper_cite = (
        "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3"
    )
    logging.info(f"\n{paper_cite}\n")
    logging.info(f"Reading database from: {DB_PATH}")
    db_name = hf.get_DB_name(DB_PATH)

    drugs = list(set(hf.get_drugs_DrugBank(DB_PATH)))

    sider_SE_dict = hf.get_side_effects_sider(
        "./DB/Data/cross_side_information_DB/SIDER/meddra_all_se.tsv"
    )
    drugBank_to_pbC, drugBank_to_pbS = hf.get_Biosnap_to_Pubchem(drugs)
    pubchem_to_stitch = hf.get_STITCH_from_Pubchem(
        set(drugBank_to_pbC.values()), set(drugBank_to_pbS.values())
    )

    drug_side_effects = []
    all_side_effects = []
    CID_PATTERN = re.compile(r"(?<=CID[ms]{1})[\d]+")
    for drug in tqdm(drugs, desc="Assigning Side Effects to drugs"):
        pubchem_id = (
            drugBank_to_pbC.get(drug, None)
            if drug in drugBank_to_pbC
            else drugBank_to_pbS.get(drug, None)
        )
        stitch = pubchem_to_stitch.get(pubchem_id, None)
        if stitch:
            stitch = CID_PATTERN.search(stitch).group()
            sider_SE = sider_SE_dict.get(stitch)
            if sider_SE:
                all_side_effects.append(sider_SE)
                drug_side_effects.append((drug, sider_SE))
            else:
                logging.debug(f"{drug} with no SE")
        else:
            logging.debug(f"{drug} not found in the SIDER database")

    drug_side_effects = dict(drug_side_effects)

    all_side_effects_positions = [
        element for sublist in all_side_effects for element in sublist
    ]
    all_side_effects_positions = sorted(list(set(all_side_effects_positions)))
    all_side_effects_positions = {
        v: k for k, v in enumerate(all_side_effects_positions)
    }

    sider_drugs_bin = []
    for drug in sorted(drugs):
        side_effects = drug_side_effects.get(drug)
        drug_se = [0] * len(all_side_effects_positions)
        if not side_effects:
            sider_drugs_bin.append(drug_se)
        else:
            for se in side_effects:
                pos = all_side_effects_positions.get(se)
                drug_se[pos] = 1
            sider_drugs_bin.append(drug_se)

    sider_bit = pd.DataFrame(
        sider_drugs_bin,
        columns=list(all_side_effects_positions.keys()),
        index=sorted(drugs),
    )
    sider_bit = pd.DataFrame(
        cosine_similarity(sider_bit), columns=sorted(drugs), index=sorted(drugs)
    )
    np.fill_diagonal(sider_bit.values, 1)

    path = hf.check_and_create_folder(db_name,model_name)
    sider_bit.to_csv(os.path.join(path, "DrugBank_drug_SIDER_SideEffect.tsv"), sep="\t")
    sider_bit.to_pickle(os.path.join(path, "DrugBank_drug_SIDER_SideEffect.pickle"))
	
def YAMANASHI_SIDER(subdataset='E', model_name = 'DDR'):
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=logging.DEBUG)

	DB_PATH = './DB/Data/Yamanashi_et_al_GoldStandard/'+subdataset+'/interactions/'+subdataset.lower()+'_admat_dgc_mat_2_line.txt'

	# sanity check for the DB
	paper_cite = 'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3'
	logging.info(f'\n{paper_cite}\n')
	logging.info(f'Reading database from: {DB_PATH}')
	db_name = hf.get_DB_name(DB_PATH)

	drugs = list(set(hf.read_and_extract_drugs_Yamanashi(DB_PATH)))

	sider_SE_dict = hf.get_side_effects_sider('./DB/Data/cross_side_information_DB/SIDER/meddra_all_se.tsv')
	kegg_SIDER_dict = hf.get_Kegg_SIDER_dict('./DB/Data/cross_side_information_DB/STITCH/dict_STITCH_KEGG.tsv')

	drug_side_effects = []
	all_side_effects = []
	for drug in drugs:
		kegg_stitch_ID = kegg_SIDER_dict.get(drug, None)
		if kegg_stitch_ID:
			sider_SE = sider_SE_dict.get(kegg_stitch_ID, None)
			if sider_SE:
				all_side_effects.append(sider_SE)
				drug_side_effects.append((drug, sider_SE))
			else:
				logging.debug(f'{drug} with no SE')	
		else:
			logging.debug(f'{drug} not found in the SIDER database')

	drug_side_effects = dict(drug_side_effects)

	all_side_effects_positions = [element for sublist in all_side_effects for element in sublist]
	all_side_effects_positions = sorted(list(set(all_side_effects_positions)))
	all_side_effects_positions = {v:k for k,v in enumerate(all_side_effects_positions)}

	sider_drugs_bin = []
	for drug in sorted(drugs):
		side_effects = drug_side_effects.get(drug)
		drug_se = [0]*len(all_side_effects_positions)
		if  not side_effects:
			sider_drugs_bin.append(drug_se)
		else:
			for se in side_effects:
				pos = all_side_effects_positions.get(se)
				drug_se[pos] = 1
			sider_drugs_bin.append(drug_se)

	sider_bit =  pd.DataFrame(sider_drugs_bin,  columns=list(all_side_effects_positions.keys()),index=sorted(drugs))

	sider_bit =  pd.DataFrame(cosine_similarity(sider_bit), columns=sorted(drugs), index=sorted(drugs))
	np.fill_diagonal(sider_bit.values, 1)

	path = hf.check_and_create_folder(db_name,model_name)
	#sider_bit.to_pickle(os.path.join(path, hf.get_yamanashi_subDB(db_name) + '_drug_SIDER_SideEffect.pickle'))
	sider_bit.to_csv(os.path.join(path, hf.get_yamanashi_subDB(db_name) + '_drug_SIDER_SideEffect.tsv'), sep='\t')
	