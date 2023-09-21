import os
import numpy as np
import pandas as pd
import logging
import argparse
import helper_functions_dtinet as hf
import json



######################################## START MAIN #########################################
#############################################################################################

def main():
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=4,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
    parser.add_argument("dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)

    args = parser.parse_args()
    # -) log info; 
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
    #######
    logging.info("============== DTIs Yamanishi ==============")
    ### log output detals
    logging.info(
        '''
        This script needs:
            - Drug Target Interactions for Davis et al
        Returns:
            - DTI tsv file --- only? 
        '''
        )

    DB_PATH = args.dbPath
    logging.debug(f'DB_PATH=={DB_PATH}')
    logging.info(f'Working in output folder for: {DB_PATH}')
    db_name = hf.get_DB_name(DB_PATH)
    hf.check_and_create_folder(db_name)
    # Create relative output path
    wdir = os.path.join('../Data', db_name)
    logging.debug(f'working directory is: {wdir}')

    db_file_path = '../../DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al_w_labels.tsv'
    davis = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Label'])
    davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
    logging.debug(davis.head(2))
    mod_davis = davis.copy()
    # we only need the positive edges for writing a edge list
    # then drop those with 0
    logging.debug('Drop null int for edge list')
    logging.debug(f'Shape before removing Label ==0 {mod_davis.shape}')
    mod_davis = mod_davis.drop(mod_davis[mod_davis.Label == 0].index)
    logging.debug(f'Shape after removing 0s {mod_davis.shape}')
    mod_davis = mod_davis.drop(columns='Label')

    ## DRUG IDENTIFIERS
    logging.debug("Changing from PubChemID to DrugBankID")
    
    #  Pubchem fmt
    mod_davis.loc[:, 'PubChemID'] = mod_davis.loc[:, 'PubChemID'].astype(int) # not float
    mod_davis.loc[:, 'PubChemID'] = mod_davis.loc[:, 'PubChemID'].astype(str) # not float
    logging.debug(mod_davis.head(2))

    logging.debug("Checking in dictionary parsing DrugBank")
    dic_cid_dbid = hf.pubchem_to_drugbankid()
    mod_davis['DrugBankID'] = mod_davis['PubChemID'].map(dic_cid_dbid)
    logging.debug(mod_davis.head(3))
    drugs_before_request = len(mod_davis.DrugBankID.unique())
    logging.debug(f'Loosing {len(mod_davis.PubChemID.unique()) - len(mod_davis.DrugBankID.unique())} drugs')

    logging.debug('Checking if we can make a comparison through SID->KEGG->DrugBank')
    list_cids_wo_DB  = mod_davis[mod_davis['DrugBankID'].isna() == True].PubChemID.unique().tolist()
    logging.debug('Creating a dictionary from CID PubChem to SID PubChem')
    cid2sid = hf.get_dic_cid2sid(list_cids_wo_DB)
    logging.debug('Creating a dictionary from KEGG to DrugBank')
    dic_kegg_db = hf.get_dict_kegg2db()
    logging.debug('Creating a dictionary from KEGG to PubChem SID')
    dic_kegg_pubchem = hf.get_dict_kegg2pubchemsid() 
    logging.debug('Trying...')
    #list_cids_wo_DB # CID
    d1 ={} 
    dic_pubchem2kegg = dict(zip(dic_kegg_pubchem.values(), dic_kegg_pubchem.keys()))
    for drug_cid in list_cids_wo_DB:
        if drug_cid in cid2sid.keys():
            for element in cid2sid.get(drug_cid):
                if element in dic_kegg_pubchem.values():
                    kegg = dic_pubchem2kegg.get(str(element))
                    drugbank_id = dic_kegg_db.get(kegg)
                    if drugbank_id:
                        d1[drug_cid] = drugbank_id

    found_d1 = list(d1.keys())
    logging.debug(f'found {len(found_d1)} drugs!')
    list_cids_wo_DB = [x for x in list_cids_wo_DB if x not in found_d1]

    logging.debug('Checking if we can make a comparison through Synonym & DrugBankID')
    logging.debug('Creating a dictionary from CID PubChem to synonyms for each compound')
    cid2syn = hf.get_dic_cid2syn(list_cids_wo_DB)
    logging.debug('Creating a dictionary from DrugBank Drug ID to Drug name')
    drugbank_dic = hf.drugname_drugbankid()
    # make sure all is lower case for true comparison
    drugbank_dic =  {k.lower(): v for k, v in drugbank_dic.items()}
    cid2syn = {k: [i.lower() for i in v] for k, v in cid2syn.items()} 

    logging.debug('Working...')
    d2 ={} 
    for drug_cid in list_cids_wo_DB:
        if drug_cid in cid2syn.keys():
            for element in cid2syn.get(drug_cid):
                if element in drugbank_dic.keys():
                    drugbank_id = drugbank_dic.get(element)
                    if drugbank_id:
                        d2[drug_cid] = drugbank_id
    found_d2 = list(d2.keys())
    logging.debug(f'found {len(found_d2)} drugs!')

    #
    logging.debug('Now we can NOW apply the maps...')
    mod_davis['DrugBankID_A'] = mod_davis['PubChemID'].map(d1)
    mod_davis['DrugBankID_B'] = mod_davis['PubChemID'].map(d2)

    mod_davis.DrugBankID.fillna(mod_davis.DrugBankID_A, inplace=True)
    mod_davis.DrugBankID.fillna(mod_davis.DrugBankID_B, inplace=True)

    drugs_after_resquest = len(mod_davis.DrugBankID.unique())
    logging.debug(f'we did our best and rescued {drugs_after_resquest-drugs_before_request} drugs')


    mod_davis = mod_davis.drop(columns=['PubChemID', 'DrugBankID_A', 'DrugBankID_B'])
    mod_davis = mod_davis.dropna()

    #### PROTEIN IDENTIFIERS
    ####### this step is not needed for Binding DB
    logging.info('We need to change from GeneName to UniprotID')
    logging.debug('Using bioMART to create a dictionary...')
    data_path_biomart = '../../DB/Data/cross_side_information_DB/bioMART/mart_export_expanded.txt'
    biomart = pd.read_csv(data_path_biomart, sep='\t', usecols=['UniProtKB Gene Name symbol', 'UniProtKB Gene Name ID'])
    biomart = biomart.dropna().drop_duplicates()
    #biomart.shape # (69394, 2)
    genename2geneid = dict(zip(biomart['UniProtKB Gene Name symbol'].tolist(), biomart['UniProtKB Gene Name ID'].tolist()))

    logging.debug('mapping...')
    mod_davis['UniprotID'] = mod_davis['GeneName'].map(genename2geneid)
    mod_davis = mod_davis.dropna().drop_duplicates()

    mod_davis = mod_davis.drop(columns='GeneName')
    
    logging.info("Now we are ready to extract information!")
    dtis_davis = mod_davis[['DrugBankID', 'UniprotID']]
    #list_of_drug_nodes = mod_davis.DrugBankID.unique().tolist()
    dict_drugid_smiles = dict(zip( mod_davis.DrugBankID.tolist(), mod_davis.SMILES.tolist() ))
    # check if all can retrieve a fp later
    dict_protein_sequence = dict(zip( mod_davis.UniprotID.tolist(), mod_davis.Sequence.tolist() ))

    logging.debug(mod_davis.head(4))
    logging.debug(f'Thera are {len(mod_davis.DrugBankID.unique())} unique drugs')
    logging.debug(f'Thera are {len(mod_davis.UniprotID.unique())} unique proteins')
    logging.info("Saving files")
    
    dtis_davis.to_csv(os.path.join(wdir, f'DTI_{DB_PATH}.tsv'), header=True,index=False ,sep="\t")

    file_path_dic_protein_seq = os.path.join(wdir, 'dic_protein_seq.json')
    with open(file_path_dic_protein_seq, 'w', encoding='utf-8') as f:
        json.dump(dict_protein_sequence, f, ensure_ascii=False, indent=4)

    file_path_dic_drug_smiles=  os.path.join(wdir, 'dic_smiles.json')
    with open(file_path_dic_drug_smiles, 'w', encoding='utf-8') as f:
        json.dump(dict_drugid_smiles, f, ensure_ascii=False, indent=4)

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-----------------------------------------------------------------------------------------
####################### END OF THE CODE ######################################################

