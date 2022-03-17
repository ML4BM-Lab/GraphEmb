import os, sys
import argparse
import logging
import numpy as np
import pandas as pd
import json
from re import search


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



######################################## START MAIN #########################################
#############################################################################################

def main():
    '''
    From HPRD:
        - protein.txt 
         - protein_dict_map.txt (or json)
        - mat_protein_protein.txt
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
    ### log output detals
    logging.info(
        '''
        This script needs:
            - XXXXX.XX
            - HUMAN_9606_idmapping.dat

        This script generates the following files:
            - genes_data.dat
            - all_protein.txt 
            - protein_dict_map.txt (or json)
            - coordinates_PPI.tsv
        
        DTINet uses from HPRD:
            - protein.txt 
            - protein_dict_map.txt
            - mat_protein_protein.txt
        '''
        )
    
    # OUTPUT DIRECTORY
    # sanity check
    DB_PATH = args.dbPath
    logging.info(f'Working in output folder for: {DB_PATH}')
    db_name = get_DB_name(DB_PATH)
    check_and_create_folder(db_name)
    # Create relative output path
    output_path = os.path.join('../Data', db_name)
    
    ## HPRD DATA FOLDER
    data_path = '../../../Data/cross_side_information_DB/HPRD'
    logging.info(f'Processing HPRD in {output_path}')

    # Read mapping data 
    logging.info(f'Reading HUMAN_9606_idmapping ...')
    with open(os.path.join(data_path,'HUMAN_9606_idmapping.dat' ), 'r') as infl:
        all_lines = infl.read().splitlines()

    all_lines = list(filter(lambda line: ('Gene_Name'  in line) or ('GeneCards' in line), all_lines))
    all_lines_spt = [i.split('\t') for i in all_lines] 
    data = pd.DataFrame(all_lines_spt, columns=['Uniprot_ID', 'type', 'GeneSymbol'])
    # test_dt.equals(data) IS TRUE ==> replace !!
    
    # get the data for processing
    # GeneName appears always before GeneCards, so its safe to do:
    translator = data.drop_duplicates(subset='GeneSymbol', keep="first") 
    # get value list for dictionary
    val_protein_id = translator['Uniprot_ID'].values.tolist()
    key_gen_name = translator['GeneSymbol'].values.tolist()
    # Create the dictionary
    dic_gen_to_protein = dict(list(zip(key_gen_name,val_protein_id)))
    # output dic
    # Write dictionary (txt file) (drugbank_ID and name)
    logging.info(f'Writing protein_dic_map.txt...')
    with open(os.path.join(output_path,'protein_dic_map.txt'), 'w') as f:
        for i in range(len(val_protein_id)):
            _ = f.write("%s:%s\n" % (val_protein_id[i],key_gen_name[i]))
    if args.json == True:
        logging.info(f'Writing protein_dic_map as json...')
        dic_protein_to_gen = dict(list(zip(val_protein_id,key_gen_name)))
        file_name_json = 'dic_protein_to_gen.json'
        with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
            json.dump(dic_protein_to_gen, f, ensure_ascii=False, indent=4)

    ## GET PROTEIN PROTEIN INTERACTIONS (coordinates)
    logging.info(f'Getting PPI coordinates...')
    hpdr_file = os.path.join(data_path, 'HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt')
    hpdr_head = ['Prot1 GeneSymbol', 'Prot1 HPRDid', 'Prot1 RefSeqid',
        'Prot2 GeneSymbol', 'Prot2 HPRDid', 'Prot2 RefSeqid',
        'Experiment type', 'Pubmed id']
    hprd = pd.read_table(hpdr_file, names=hpdr_head, usecols=['Prot1 GeneSymbol', 'Prot1 HPRDid','Prot2 GeneSymbol', 'Prot2 HPRDid'])
    # map from gen to protein
    hprd['Prot1 UniprotKB'] = hprd['Prot1 GeneSymbol'].map(dic_gen_to_protein)
    hprd['Prot2 UniprotKB'] = hprd['Prot2 GeneSymbol'].map(dic_gen_to_protein)
    # drop HPDR (dont actually need tat)
    hprd_processed =  hprd.drop(columns=['Prot1 HPRDid','Prot2 HPRDid'])
    # first drop rows that contain no Gene Symbol ('-')
    hprd_processed = hprd_processed.drop(index=hprd_processed[hprd_processed['Prot1 GeneSymbol']=='-'].index, axis=1)
    hprd_processed = hprd_processed.drop(index=hprd_processed[hprd_processed['Prot2 GeneSymbol']=='-'].index, axis=1)
    # drop nan
    hprd_processed = hprd_processed.dropna()
    hprd_processed[hprd_processed['Prot1 UniprotKB'].isna()]
    hprd_processed[hprd_processed['Prot2 UniprotKB'].isna()]
    # this can be saved as: (uncoment)
    # hprd_processed.to_csv(os.path.join(output_path,'hprd_processed_with_uniprot.tsv'), sep="\t")
    # Get a DataFrame:
    df_PPI = hprd_processed.drop(['Prot1 GeneSymbol', 'Prot2 GeneSymbol'], axis=1)
    df_PPI.columns = ['Prot1', 'Prot2']
    df_PPI = df_PPI.drop_duplicates()
    # save coordinate files
    logging.info(f'Writting PPI coordinate file...')
    df_PPI.to_csv(os.path.join(output_path, 'coordinates_PPI.tsv'), index=False, header=True, sep='\t')

    ########### 
    ## GET PROTEIN NODES
    s1 = hprd_processed['Prot1 UniprotKB'].values.tolist()
    s2 = hprd_processed['Prot2 UniprotKB'].values.tolist()
    protein_list = list(set(s1 + s2))
    logging.info(f'Writing all_protein.txt file...')
    np.savetxt(os.path.join(output_path, 'all_protein.txt'), protein_list, newline='\n', fmt='%s')
    
    # remove genes_data.dat before finishing
    # os.system(f'rm {output_path}/genes_data.dat')
    logging.info(f'Done with HPRD!')




#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
