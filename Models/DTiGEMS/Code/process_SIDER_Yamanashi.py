import helper_functions_dtigems as hf
from sklearn.metrics.pairwise import cosine_similarity
import logging
import argparse
import pandas as pd
import numpy as np
import os

######################################## START MAIN #########################################
#############################################################################################
def main():
    # get the parameters from the user
    parser = argparse.ArgumentParser()
    parser.add_argument("dbPath", help="Path to the database", type=str)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbosity",
        action="count",
        default=3,
        help="Verbosity (between 1-4 occurrences with more leading to more "
        "verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
        "DEBUG=4",
    )
    args = parser.parse_args()
    # define logging level
    log_levels = {
        0: logging.CRITICAL,
        1: logging.ERROR,
        2: logging.WARN,
        3: logging.INFO,
        4: logging.DEBUG,
    }

    # set the logging info
    level = log_levels[args.verbosity]
    fmt = "[%(levelname)s] %(message)s"
    logging.basicConfig(format=fmt, level=level)

    # sanity check for the DB
    paper_cite = (
        "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3"
    )
    logging.info(f"\n{paper_cite}\n")
    DB_PATH = args.dbPath
    # DB_PATH =  './../../DB/Data/Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt'
    logging.info(f"Reading database from: {DB_PATH}")
    db_name = hf.get_DB_name(DB_PATH)

    drugs = list(set(hf.read_and_extract_drugs_Yamanashi(DB_PATH)))

    sider_SE_dict = hf.get_side_effects_sider(
        "./../../DB/Data/cross_side_information_DB/SIDER/meddra_all_se.tsv"
    )
    kegg_SIDER_dict = hf.get_Kegg_SIDER_dict(
        "./../../DB/Data/cross_side_information_DB/STITCH/dict_STITCH_KEGG.tsv"
    )

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

    path = hf.check_and_create_folder(db_name)
    sider_bit.to_csv(
        os.path.join(
            path, hf.get_yamanashi_subDB(db_name) + "_Drug_SIDER_SideEffect.tsv"
        ),
        sep="\t",
    )
    sider_bit.to_pickle(
        os.path.join(
            path, hf.get_yamanashi_subDB(db_name) + "_Drug_SIDER_SideEffect.pickle"
        )
    )


#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
    main()

#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
