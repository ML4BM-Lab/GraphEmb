import os
import numpy as np
import pandas as pd
import logging
import argparse
import helper_functions_dtinet as hf



######################################## START MAIN #########################################
#############################################################################################

def main():
    parser = argparse.ArgumentParser() 
    parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=3,
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
            - Drug Target Interactions for each type (NR, IC, GPCR, E)
        Returns:
            - DTI tsv file
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
    colnames_data = ['Kegg_ID', 'Gene']
    data_path = '../../DB/Data/Yamanashi_et_al_GoldStandard/' + DB_PATH.upper()+'/interactions/' + DB_PATH.lower() +'_admat_dgc_mat_2_line.txt'
    df = pd.read_csv(data_path, header = None, names = colnames_data, index_col=False, sep='\t')
    df['Gene'] = df['Gene'].map(lambda x: x[:3] + ":" + x[3:])
    logging.debug(f'{(df.head(2))}')

    #
    logging.info('Loading hsa2uni dictionary...')
    hsa2uni = hf.get_dict_hsa2uni()
    # repeated entries, but diff uniprots have the same sequence (so it's safe)
    df['Uniprot'] = df['Gene'].map(hsa2uni)
    logging.debug(f'{(df.head(2))}')
    # Change from KEGG Drug ID to DrugBank  
    logging.debug("Loading kegg2db dictionary...")
    kegg2db = hf.get_dict_kegg2db()
    df['DrugBank_ID'] = df['Kegg_ID'].map(kegg2db)
    # Process df before saving
    logging.info("Processing & saving DTIs to file...")
    df = df.drop(columns=['Gene', 'Kegg_ID']) #
    df.columns = ['Protein', 'DrugBank_ID']
    DTI = df[['DrugBank_ID', 'Protein']] 
    # DTI.isna().sum()
    DTI = DTI.dropna()
    DTI.drop_duplicates()
    logging.debug(f'{(df.head(2))}')
    logging.info(f'    Matrix Shape: {DTI.shape}')
    logging.info(f'    Unique drugs: {len(DTI.DrugBank_ID.unique())}')
    logging.info(f'    Unique targets: {len(DTI.Protein.unique())}')
    #DTI
    DTI.to_csv(os.path.join(wdir, f'DTI_{DB_PATH}.tsv'), header=True,index=False ,sep="\t")


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-----------------------------------------------------------------------------------------
####################### END OF THE CODE ######################################################


