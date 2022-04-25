import os
import numpy as np
import pandas as pd
import logging
import argparse
import helper_functions_dtinet as hf
import json
from tqdm import tqdm

########### Modify Davis for BindingDB --------->>>>> 

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

    db_file_path = '../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
    bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Y'])
    bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID', 'Target Sequence': 'Sequence'}, axis=1)


    logging.debug("Changing from PubChemID to DrugBankID")

    mod_bind = bindingdb.copy()

    threshold = 30

    mod_bind['Label'] = [1 if x < threshold else 0 for x in mod_bind['Y']]
    mod_bind = mod_bind.drop(columns='Y')

    logging.debug('Drop null int for edge list')
    logging.debug(f'Shape before removing Label ==0 {mod_bind.shape}')
    mod_bind = mod_bind.drop(mod_bind[mod_bind.Label == 0].index)
    logging.debug(f'Shape after removing 0s {mod_bind.shape}')
    mod_bind = mod_bind.drop(columns='Label')

    # mod_bind = mod_bind.drop(columns=['SMILES', 'GeneName', 'Sequence'])
    ##### work with drugs
    mod_bind.loc[:, 'PubChemID'] = mod_bind.loc[:, 'PubChemID'].astype(int) # not float
    mod_bind.loc[:, 'PubChemID'] = mod_bind.loc[:, 'PubChemID'].astype(str) # not float
    logging.debug(mod_bind.head(2))


    logging.debug("Checking in dictionary parsing DrugBank")
    dic_cid_dbid = hf.pubchem_to_drugbankid()
    mod_bind['DrugBankID'] = mod_bind['PubChemID'].map(dic_cid_dbid)
    mod_bind.head()
    drugs_before_request = len(mod_bind.DrugBankID.unique())
    logging.debug(f'Loosing {len(mod_bind.PubChemID.unique()) - len(mod_bind.DrugBankID.unique())} drugs')

    logging.debug('Checking if we can make a comparison through SID->KEGG->DrugBank')
    list_cids_wo_DB  = mod_bind[mod_bind['DrugBankID'].isna() == True].PubChemID.unique().tolist()
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
    for drug_cid in tqdm(list_cids_wo_DB, desc='Checking SID->KEGG->DrugBank', position=0, leave=True):
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
    for drug_cid in tqdm(list_cids_wo_DB, desc='Checking CID->Name->DrugBank', position=0, leave=True):
        if drug_cid in cid2syn.keys():
            for element in cid2syn.get(drug_cid):
                if element in drugbank_dic.keys():
                    drugbank_id = drugbank_dic.get(element)
                    if drugbank_id:
                        d2[drug_cid] = drugbank_id

    found_d2 = list(d2.keys())
    logging.debug(f'found {len(found_d2)} drugs!')

    ###### working here --------->>> ** * 
    logging.debug('Now we can NOW apply the maps...')
    mod_bind['DrugBankID_A'] = mod_bind['PubChemID'].map(d1)
    mod_bind['DrugBankID_B'] = mod_bind['PubChemID'].map(d2)

    mod_bind.DrugBankID.fillna(mod_bind.DrugBankID_A, inplace=True)
    mod_bind.DrugBankID.fillna(mod_bind.DrugBankID_B, inplace=True)

    drugs_after_resquest = len(mod_bind.DrugBankID.unique())
    logging.debug(f'we did our best and rescued {drugs_after_resquest-drugs_before_request} drugs')


    mod_bind = mod_bind.drop(columns=['PubChemID', 'DrugBankID_A', 'DrugBankID_B'])
    mod_bind = mod_bind.dropna()



    logging.info("Now we are ready to extract information!")
    dtis_davis = mod_bind[['DrugBankID', 'UniprotID']]
    #list_of_drug_nodes = mod_bind.DrugBankID.unique().tolist()
    dict_drugid_smiles = dict(zip( mod_bind.DrugBankID.tolist(), mod_bind.SMILES.tolist() ))
    # check if all can retrieve a fp later
    dict_protein_sequence = dict(zip( mod_bind.UniprotID.tolist(), mod_bind.Sequence.tolist() ))

    logging.debug(mod_bind.head(4))
    logging.debug(f'Thera are {len(mod_bind.DrugBankID.unique())} unique drugs')
    logging.debug(f'Thera are {len(mod_bind.UniprotID.unique())} unique proteins')
    logging.info("Saving files...")
    
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

