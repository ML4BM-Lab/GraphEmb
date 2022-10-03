import argparse
import os
import helper_functions_DTI2Vec as hf
import logging
import subprocess as sp
from tqdm import tqdm


######################################## START MAIN #########################################
#############################################################################################
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dbPath",
        help="Path to the database interaction lits",
        default="./../../DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv",
        type=str,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbosity",
        default="3",
        help="Verbosity (between 1-4 occurrences with more leading to more "
        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
        "DEBUG=4",
    )
    args = parser.parse_args()
    log_levels = {
        "0": logging.CRITICAL,
        "1": logging.ERROR,
        "2": logging.WARN,
        "3": logging.INFO,
        "4": logging.DEBUG,
    }
    # set the logging info
    level = log_levels[args.verbosity]
    fmt = "[%(levelname)s] %(message)s"
    logging.basicConfig(format=fmt, level=level)
    # sanity check for the DB
    # DB_PATH = './../../DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
    DB_PATH = args.dbPath
    logging.info(f"Reading database from: {DB_PATH}")
    db_name = hf.get_DB_name(DB_PATH)
    for K_neigh in [2, 5, 10, 20]:
        logging.info(f"Running for K={K_neigh}")
        dtis, node_index_dict = hf.read_dtis_BIOSNAP(DB_PATH)
        new_drug_edges = hf.get_k_neighbors(f'./../Data/{db_name}/Drug_simmatrix/simmatrix/TFIDF__media_scratch_ssd_tmp_BIOSNAP_simmat_dc.txt', top_k=K_neigh)
        dtis.extend(hf.create_edges(new_drug_edges))
        new_prot_edges = hf.get_k_neighbors(f'./../Data/{db_name}/Proteins_SmithWaterman_scores_MinMax.tsv', top_k=K_neigh)
        dtis.extend(hf.create_edges(new_prot_edges, replace_dots=True))
        admat = hf.get_admat_from_dti(dtis)
        file_path = os.path.join('./../Data',db_name , db_name + '_' +str(K_neigh) )
        logging.debug(f"Output files at : {file_path}")
        admat.to_csv(file_path +'_admat.tsv', sep='\t')
        hf.write_edges(dtis, file_path +'_dti.tsv')
        # encode the DTIs
        logging.info("Encoding the DTIs")
        encoded_dtis = [
            (node_index_dict.get(node1), node_index_dict.get(node2))
            for node1, node2 in dtis
        ] 
        TMP_PATH = hf.create_remove_tmp_folder(os.path.join("/tmp/N2V", db_name))

        dti_coded_path = hf.write_dtis(encoded_dtis, TMP_PATH)
        dict_path = os.path.splitext(dti_coded_path)[0] + "_node_index_dict.txt"
        hf.write_dict(node_index_dict, dict_path)

        PATH_N2V = "/home/margaret/data/gserranos/SNAP/examples/node2vec/node2vec"
        # The node_2_vec implementation needs the names of the nodes as numeric elements, so we need to build dictionaries
        dimensions = [64, 128, 256, 720, 1440]
        for cv in range(10):
            cv = cv + 1
            for dim in tqdm(dimensions, desc="Creating the embeddings"):
                embeddings_file = (
                    f"./../Data/{db_name}/BIOSNAP_Nodes_Embedding_{str(dim)}_{cv}.emb"
                )
                cmd = f"{PATH_N2V} -i:{dti_coded_path} -o:{embeddings_file} -d:{dim} -v"
                out_1 = sp.check_call(cmd, shell=True, stdout=sp.DEVNULL)
                if out_1 == 0:
                    hf.write_embedddings_with_name(embeddings_file, node_index_dict)
                else:
                    logging.error(
                        f"Something went wrong with the node2vec implementation. Check the output of the command: {cmd}"
                    )


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
