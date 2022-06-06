import argparse
import os
import helper_functions_DTI2Vec as hf
import logging
import subprocess as sp
from tqdm import tqdm


######################################## START MAIN #########################################
#############################################################################################




def get_k_neighbors(path, top_k=5):
    # './../Data/Yamanashi_et_al_GoldStandard/NR/Drugs_SIMCOMP_scores.tsv'
    with open(path, 'r') as f:
        drug_drug_sim = f.readlines()
    #
    new_edges=[]
    drug_drug_sim = [distance_vec.strip().split('\t') for distance_vec in drug_drug_sim]
    drug_drug_sim_index = drug_drug_sim.pop(0)
    drug_drug_sim_metric = [drug_sim[1:] for drug_sim in drug_drug_sim]
    drug_drug_sim_metric = [list(map(float, line)) for line in drug_drug_sim_metric]
    #
    for drug in range(len(drug_drug_sim_index)):
        neighbours = sorted(drug_drug_sim_metric[drug], reverse=True)
        # +1 to keep trully the top_k, not himself
        neighbours = neighbours[:(top_k +1)]
        neighbour_names = []
        for neigh in neighbours[1:]:
            neigh_idx = drug_drug_sim_metric[drug].index(neigh)
            neigh_name = drug_drug_sim_index[neigh_idx]
            if neigh_name not in neighbour_names:
                neighbour_names.append(neigh_name)
            else:
                neigh_name = drug_drug_sim_index[drug_drug_sim_metric[drug].index(neigh, neigh_idx+1)]
                neighbour_names.append(neigh_name)
        new_edges.append((drug_drug_sim_index[drug], neighbour_names))
    return new_edges

def create_edges(edge_list):
    new_edges = []
    for orig, dest_list in edge_list:
        new_edges.extend([[orig, dest] for dest in dest_list])
    return new_edges

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dbPath",
        help="Path to the database interaction lits",
        default="/home/margaret/data/jfuente/DTI/Data/Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt",
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
    fmt = "[%(levelname)s] %(message)s]"
    logging.basicConfig(format=fmt, level=level)

    # sanity check for the DB
    # DB_PATH = './../../DB/Data/Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt'
    DB_PATH = args.dbPath
    logging.info(f"Reading database from: {DB_PATH}")
    db_name = hf.get_DB_name(DB_PATH)
    dtis, node_index_dict = hf.read_dtis_Yamanashi(DB_PATH)

    # TODO: check index names, not working with the current version of the code
    #data augmentation
    new_drug_edges = get_k_neighbors('./../Data/Yamanashi_et_al_GoldStandard/NR/Drugs_SIMCOMP_scores.tsv', top_k=5)
    dtis.extend(create_edges(new_drug_edges))
    new_prot_edges = get_k_neighbors('./../Data/Yamanashi_et_al_GoldStandard/NR/Proteins_SmithWaterman_scores_MinMax.tsv', top_k=5)
    dtis.extend(create_edges(new_prot_edges))
    # encode the DTIs
    logging.info("Encoding the DTIs")
    encoded_dtis = [
        (node_index_dict.get(node1), node_index_dict.get(node2))
        for node1, node2 in dtis
    ]
    TMP_PATH = hf.create_remove_tmp_folder(os.path.join('/media/scratch_ssd/tmp/', db_name))

    dti_coded_path = hf.write_dtis(encoded_dtis, TMP_PATH)
    dict_path = os.path.splitext(dti_coded_path)[0] + "_node_index_dict.txt"
    hf.write_dict(
        node_index_dict, os.path.splitext(dti_coded_path)[0] + "_node_index_dict.txt"
    )

    PATH_N2V = "/home/margaret/data/gserranos/SNAP/examples/node2vec/node2vec"
    # The node_2_vec implementation needs the names of the nodes as numeric elements, so we need to build dictionaries
    dimensions = [64, 128, 256, 720, 1440]
    for cv in range(10):
        cv = cv + 1
        subname = hf.get_yamanashi_subDB(db_name)
        for dim in tqdm(dimensions, desc="Creating the embeddings"):
            embeddings_file = (
                f"./../Data/{db_name}/{subname}_Nodes_Embedding_{str(dim)}_{cv}.emb"
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
