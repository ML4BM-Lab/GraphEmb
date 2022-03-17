import os, sys
import xml.etree.ElementTree as ET
from tqdm import tqdm
import argparse
import logging
from tqdm import tqdm
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



#################### GET COORDINATES of DTI INTERACTIONS #################################### ----> esta funcion está mal porque non ten en conta enzymes, transport & carrriers
############################################################################################# ---->
def get_DTI_coordinates(drug_entry, coordinate_list):
	drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
	target = 'None'
	for tgt in  drug_entry.findall('.//{http://www.drugbank.ca}targets')[0]: # iterate in targets
		for i in tgt.findall('.//{http://www.drugbank.ca}polypeptide'): # for searching external-identifiers
			ext = i.findall('.//{http://www.drugbank.ca}external-identifiers')[0]
			for id in list(ext):
				if list(id)[0].text == "UniProtKB": # busca donde esta uniprot
					target = list(id)[1].text
					if target != 'None':
						coordinate = (drugbank_ID, target)
						coordinate_list.append(coordinate)
					else:
						continue 
    # enzymes
    # transporters
    # carrier
	return coordinate_list


######################################## START MAIN #########################################
#############################################################################################

def main():
    '''
    DTIs from DrugBank
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
    ### log output detals
    logging.info(
        '''
        This script needs:
            - DrugBank xml file

        This script generates:
            - coordinates_DTI.tsv
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

    ## program
    logging.info(f'Reading DrugBank xml file')
    tree = ET.parse('../../../Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    root = tree.getroot()
    # execute the function
    logging.info(f'Retrieving DTI coordinates')
    coordinate_list = []
    for drug_entry in tqdm(root):
        coordinate_list = get_DTI_coordinates(drug_entry, coordinate_list)
    # save coordinates to a tsv
    with open(os.path.join(output_path, 'coordinates_DTI.tsv'), 'w') as f:
        f.write("DrugBank_ID\tProtein_ID\n")
        for item in coordinate_list:
            _ = f.write("%s\t%s\n" % item)

#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################

