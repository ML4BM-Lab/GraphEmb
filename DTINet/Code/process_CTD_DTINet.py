import os, sys
import argparse
import logging
import numpy as np
import pandas as pd
import json
from tqdm import tqdm
from re import search
import xml.etree.ElementTree as ET


def get_DB_name(path):
	"""
	This function returns the name of the DB.
	"""
	DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank', 'E', 'GPCR', 'IC', 'NR']
	for db in DB_NAMES:
		if search(db, path):
			logging.info(f'Database: {db}')
			if db in ['E', 'GPCR', 'IC', 'NR']:
				db = os.path.join('Yamanashi_et_al_GoldStandard', db)
				return db
			else:
				return db
	logging.error(f'Database: {db} not found')
	sys.exit('Please provide a valid database')


def check_and_create_folder(db_name):
	if not os.path.exists(os.path.join('../Data', db_name)):
		os.mkdir(os.path.join('../Data', db_name))


def drugbankid_to_pubchem():
    tree = ET.parse('../../../Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    root = tree.getroot()
    dbids = []
    pubchemids = []
    for drug_entry in tqdm(root):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        pubchem_id = np.nan # to not repeat id if it does not appear later
        for props in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
            for prop in props:
                if(prop.text == 'PubChem Compound'): 
                    pubchem_id = props[1].text
                break # romper una vez que encuentre el pubchem 
        dbids.append(drugbank_ID)
        pubchemids.append(pubchem_id)

    dic_cid_dbid = dict((zip(pubchemids, dbids)))
    # first element is nan, delete
    dic_cid_dbid[list(dic_cid_dbid.keys())[0]]
    dic_cid_dbid.pop(list(dic_cid_dbid.keys())[0])
    return dic_cid_dbid


######################################## START MAIN #########################################
#############################################################################################


def main():
    '''
    CTD
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
                    help="Verbosity (between 1-4 occurrences with more leading to more "
                        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
                        "DEBUG=4")
    parser.add_argument('-json', help="If selected, outputs a dictionary in json", action="store_true")
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
    ###  output details
    logging.info(
        '''
        This script needs:
            - CTD Database (CTD_diseases.csv, CTD_genes.csv, CTD_genes_diseases.csv, CTD_chemicals_diseases.csv)
            - all_protein.txt
            - dic_drugnames_cid.json (from XXXXX.py)
        
        This script generates the following files:
            - disease.txt
            - disease_dict_map.txt (or json)
            - coordinates_protein_disease.tsv
            - XXXXXX 
        
        DTINet uses from CTD the processed files:
            - mat_protein_disease.txt
            - mat_drug_disease.txt
        '''
        )

    # OUTPUT DIRECTORY
    DB_PATH = args.dbPath
    db_name = get_DB_name(DB_PATH)
    check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)
    ## CTD DATA FOLDER
    data_path = '../../../Data/cross_side_information_DB/CTD'
    # info 
    logging.info(f'Processing CTD in {output_path}')

    #### DISEASE NODES
    logging.info(f'Using: CTD_diseases.csv')
    header_dis_dic = ['DiseaseName', 'DiseaseID', 'AltDiseaseIDs', 'Definition', 'ParentIDs', 
        'TreeNumbers', 'ParentTreeNumbers', 'Synonyms', 'SlimMappings']

    dis_dict = pd.read_csv(os.path.join(data_path, 'CTD_diseases.csv'), 
                            index_col=False, 
                            names=header_dis_dic, 
                            comment='#', 
                            usecols=['DiseaseName', 'DiseaseID'])
    dis_dict.DiseaseName = dis_dict.iloc[:, :1].applymap(lambda s: s.lower() if type(s) == str else s)

    # In case we only want the number withouth MESH: or OMIM:
    # dis_dict.DiseaseID = dis_dict.iloc[:, 1:2].applymap(lambda s: s[5:] if type(s) == str else s)

    # export to files:::: diseases.txt is a list of diseases
    logging.info(f'Writting disease.txt')
    np.savetxt(os.path.join(output_path,'disease.txt') , dis_dict.DiseaseName.values, newline='\n', fmt='%s')
    # export to files:::: disease_dict_map.txt is a disctionary of diseases and identifiers
    logging.info(f'Writting disease_dict_map.txt')
    np.savetxt(os.path.join(output_path, 'disease_dict_map.txt'), dis_dict.DiseaseName.values + ',' + dis_dict.DiseaseID.values, newline='\n', fmt='%s')
    #### IF JSON:
    if args.json == True:
        logging.info(f'Writing disease_dict_map as json...')
        dic_disease_to_disID = dict(list(zip(dis_dict.DiseaseName.values, dis_dict.DiseaseID.values)))
        file_name_json = 'dic_disease_to_disID.json'
        with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
            json.dump(dic_disease_to_disID, f, ensure_ascii=False, indent=4)

    ####################### PROTEIN DISEASE ASSOCIATIONS #######################
    logging.info('Reading CTD_genes.csv to get a dict')
    header_gene_voc = ['GeneSymbol','GeneName','GeneID','AltGeneIDs','Synonyms','BioGRIDIDs','PharmGKBIDs','UniProtIDs']
    gen_voc = pd.read_csv(os.path.join(data_path, 'CTD_genes.csv'),
                            index_col=False, names=header_gene_voc, comment='#', 
                            usecols=['GeneSymbol','UniProtIDs'])

    gen_voc = gen_voc.dropna()
    gen_voc.shape 
    # los uniprots se pueden separar por '|' y ver si estan en la lista set proteis.txt
    gen_voc.UniProtIDs = gen_voc.UniProtIDs.str.split("|")
    gen_voc = gen_voc.explode('UniProtIDs') # as different entry

    # now filter those that are not in the human protein set from HPRD (all_protein.txt) 
    logging.info('Reading all_protein.txt for filtering human proteins')
    protein_file = os.path.join(output_path, 'all_protein.txt')
    with open(protein_file, 'r') as fl:
        hum_prots = fl.read().splitlines()

    ## filter and keep only those rows that have a HUMAN UniProt ID !!
    gen_voc = gen_voc.astype({"UniProtIDs": str})
    gen_voc['comp'] = gen_voc['UniProtIDs'].isin(hum_prots)
    gen_voc[gen_voc["comp"]==True].count()  # da 9061 uf
    gene_and_uniprot = gen_voc[gen_voc["comp"]==True][['GeneSymbol','UniProtIDs']]

    # create a dictinonary for mappig later to uniprot!
    # genesymbol == keys ; uniprot ids == values
    keys_gensymb = gene_and_uniprot['GeneSymbol'].values.tolist()
    values_uniprot = gene_and_uniprot['UniProtIDs'].values.tolist()
    dic_gen_to_protein_CTD = dict(list(zip(keys_gensymb,values_uniprot)))

    #################
    ## genes diseases # use gene symbol as principal key (sq)
    gen_dis_head = ['GeneSymbol', 'GeneID', 'DiseaseName', 'DiseaseID',
                    'DirectEvidence', 'InferenceChemicalName', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

    logging.info('Reading CTD_genes_diseases.csv')
    gen_dis = pd.read_csv(os.path.join(data_path, 'CTD_genes_diseases.csv'),
                            index_col=False, names=gen_dis_head, comment='#', usecols= ['GeneSymbol', 'DiseaseID'])

    gen_dis.isnull().values.any()
    gen_dis = gen_dis.dropna()
    # gen_dis is a dataframe that relates each gene with a DiseaseID
    # now use the dictionary created above to map
    gen_dis['UniprotID'] = gen_dis['GeneSymbol'].map(dic_gen_to_protein_CTD)
    gen_dis = gen_dis.dropna()
    protein_disease = gen_dis[['UniprotID', 'DiseaseID']]
    protein_disease = protein_disease.drop_duplicates()
    logging.info(f'shape of coordinate file is: {protein_disease.shape}')
    logging.info('Writing coordinates_protein_disease.tsv...')
    protein_disease.to_csv(os.path.join(output_path, 'coordinates_protein_disease.tsv'), index=False, header=True, sep="\t")


    ####################### DRUG - DISEASE ASSOC ####################### ---------------------------------------->>> ******
    logging.info('Reading CTD_chemicals_diseases.csv')
    h_chem = ['ChemicalName', 'ChemicalID' ,'CasRN', 'DiseaseName', 'DiseaseID',
            'DirectEvidence', 'InferenceGeneSymbol', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

    chem_dis = pd.read_csv(os.path.join(data_path, 'CTD_chemicals_diseases.csv'), 
                                index_col=False, 
                                names=h_chem, 
                                comment='#',
                                usecols=['ChemicalName', 'DiseaseID', 'DiseaseID'])

    ## remove NaN
    chem_dis = chem_dis.dropna()

    # we need the DrugBankID, as is the identification that we use in DTINet
    # Chemical ID is a MESH ID that appears in Pubchem, buy as synonym
    # we cannot ask directly for a synonym to retrieve the CID
    # then change from ChemicalName to PubchemCID using pubchempy

    logging.info('Reading dic_drugnames_cid.json ...')
    PATH_dic_drugnames_cid = '../../../Data/cross_side_information_DB/dic_drugnames_cid.json'
    # read the dictionary (slow; just read but make the script above available) 
    # This fule can be obtained running: get_dic_names_to_CID.py
    with open(PATH_dic_drugnames_cid, 'r') as f:
        dic_drugnames_cid = json.load(f)

    # Then, we need to go back to DrugBank
    # and make a relation between PubChem and DrugBank
    ################## esto esta repetido en SIDER ######
    logging.info('Reading DrugBank xml file...')
    dic_cid_dbid = drugbankid_to_pubchem()

    ## apply two maps
    chem_dis['PubChemID'] = chem_dis['ChemicalName'].map(dic_drugnames_cid)
    chem_dis.shape
    chem_dis = chem_dis.dropna()
    chem_dis.shape # loosing 1158846 entries here
    chem_dis['PubChemID'] = chem_dis['PubChemID'].astype(int).astype(str)
    chem_dis['DrugBankID'] = chem_dis['PubChemID'].apply(lambda x: dic_cid_dbid[x] if x in dic_cid_dbid else np.nan)
    chem_dis = chem_dis.dropna()
    # creation of a coordinate dataframe:
    coordinates_drug_dis = chem_dis[['DrugBankID', 'DiseaseID']]
    coordinates_drug_dis = coordinates_drug_dis.drop_duplicates()
    logging.info(f'shape of coordinate file is: {coordinates_drug_dis.shape}')
    logging.info('Writing coordinate file drug-disease...')
    coordinates_drug_dis.to_csv(os.path.join(output_path, 'coordinates_drug_disease.tsv'), index=False, header=True, sep="\t")
    # in preproc remove coordinates before creating matrix!




#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################

