import os, sys
import argparse
import logging
import numpy as np
import pandas as pd
import json
from tqdm import tqdm
from re import search
import xml.etree.ElementTree as ET
import time
import pubchempy as pcp
import multiprocessing as mp
import helper_functions_dtinet as hf
#from parallelbar import progress_map

######################################## START MAIN #########################################
#############################################################################################


def main():
    '''
    CTD
    '''
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=4,
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
    logging.info("============== CTD ==============")
    # Output details
    logging.info(
        '''
        This script needs:
            - CTD Database (CTD_diseases.csv, CTD_genes.csv, CTD_genes_diseases.csv, CTD_chemicals_diseases.csv)
            - all_protein.txt
        
        This script generates the following files:
            - disease_dict_map.txt (or json)
            - edgelist_protein_disease.tsv
            - edgelist_drug_disease.tsv
        
        DTINet uses from CTD the processed files:
            - mat_protein_disease.txt
            - mat_drug_disease.txt
        '''
        )

    # OUTPUT DIRECTORY
    DB_PATH = args.dbPath
    db_name = hf.get_DB_name(DB_PATH)
    hf.check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)
    ## CTD DATA FOLDER
    data_path = '../../DB/Data/cross_side_information_DB/CTD'
    # remove info stat from pubchempy
    logger = logging.getLogger('pubchempy') #
    logger.setLevel(level = logging.DEBUG)
    logger.disabled = True
    # info 
    logging.info(f'Processing CTD in {output_path}')

    #### DISEASE NODES
    logging.debug(f'Using: CTD_diseases.csv')
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
    #np.savetxt(os.path.join(output_path,'disease.txt') , dis_dict.DiseaseName.values, newline='\n', fmt='%s')
    # export to files:::: disease_dict_map.txt is a disctionary of diseases and identifiers
    logging.debug(f'Writting disease_dict_map.txt')
    np.savetxt(os.path.join(output_path, 'disease_dict_map.txt'), dis_dict.DiseaseName.values + ',' + dis_dict.DiseaseID.values, newline='\n', fmt='%s')
    #### IF JSON:
    if args.json == True:
        logging.debug(f'Writing disease_dict_map as json...')
        dic_disease_to_disID = dict(list(zip(dis_dict.DiseaseName.values, dis_dict.DiseaseID.values)))
        file_name_json = 'dic_disease_to_disID.json'
        with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
            json.dump(dic_disease_to_disID, f, ensure_ascii=False, indent=4)

    ####################### PROTEIN DISEASE ASSOCIATIONS #######################
    logging.debug('Reading CTD_genes.csv to get a dict')
    header_gene_voc = ['GeneSymbol','GeneName','GeneID','AltGeneIDs','Synonyms','BioGRIDIDs','PharmGKBIDs','UniProtIDs']
    gen_voc = pd.read_csv(os.path.join(data_path, 'CTD_genes.csv'),
                            index_col=False, names=header_gene_voc, comment='#', 
                            usecols=['GeneSymbol','UniProtIDs'])

    gen_voc = gen_voc.dropna()
    gen_voc.shape 
    gen_voc.UniProtIDs = gen_voc.UniProtIDs.str.split("|")
    gen_voc = gen_voc.explode('UniProtIDs') # as different entry

    # now filter those that are not in the human protein set from HPRD (all_protein.txt) 
    logging.debug('Reading all_protein.txt for filtering human proteins...')
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
    logging.info(f'    Getting protein-disease coordinates...')
    ## genes diseases # use gene symbol as principal key (sq)
    gen_dis_head = ['GeneSymbol', 'GeneID', 'DiseaseName', 'DiseaseID',
                    'DirectEvidence', 'InferenceChemicalName', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

    logging.debug('Reading CTD_genes_diseases.csv ...')
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
    logging.info(f'    shape of coordinate file is: {protein_disease.shape}')
    logging.debug('Writing coordinates_protein_disease.tsv...')
    protein_disease.to_csv(os.path.join(output_path, 'edgelist_protein_disease.tsv'), index=False, header=True, sep="\t")


    ## DRUG - DISEASE ASSOC 
    logging.info(f'    Getting drug-disease coordinates...')
    logging.debug('Reading CTD_chemicals_diseases.csv ...')
    h_chem = ['ChemicalName', 'ChemicalID' ,'CasRN', 'DiseaseName', 'DiseaseID',
            'DirectEvidence', 'InferenceGeneSymbol', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

    chem_dis = pd.read_csv(os.path.join(data_path, 'CTD_chemicals_diseases.csv'), 
                                index_col=False, 
                                names=h_chem, 
                                comment='#',
                                usecols=['ChemicalName', 'DiseaseID', 'DiseaseID'])

    ## remove NaN
    chem_dis = chem_dis.dropna()
    #
    logging.debug('Reading DrugBank xml file...')
    drugbank_dic = hf.drugname_drugbankid()

    chem_dis['DrugBankID'] = chem_dis['ChemicalName'].map(drugbank_dic)
    chem_dis = chem_dis.dropna()
    # creation an edge list as dataframe:
    coordinates_drug_dis = chem_dis[['DrugBankID', 'DiseaseID']]
    coordinates_drug_dis = coordinates_drug_dis.drop_duplicates()
    logging.info(f'     shape of coordinate file is: {coordinates_drug_dis.shape}')
    logging.debug('Writing edge list file drug-disease...')
    coordinates_drug_dis.to_csv(os.path.join(output_path, 'edgelist_drug_disease.tsv'), index=False, header=True, sep="\t")


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-----------------------------------------------------------------------------------------
####################### END OF THE CODE ######################################################

