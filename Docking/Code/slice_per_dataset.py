import os
import pandas as pd
import numpy as np
import helper_functions as hf
import logging


logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)

## full matrices in pkl format
PKLS_PATH = '../Data/pkls'

# For Drugs
drugs_tani = pd.read_pickle(os.path.join(PKLS_PATH, 'drugs_tani_full.pkl'))
drugs_tani.index = drugs_tani.index.astype(str)
drugs_tani.columns = drugs_tani.columns.astype(str)
# Annotations
drugs_annot_fam = pd.read_pickle(os.path.join(PKLS_PATH, 'df_annot_fam.pkl'))
drugs_annot_fam = drugs_annot_fam.astype(str)

drugs_annot_moldesc = pd.read_pickle(os.path.join(PKLS_PATH, 'annot_moldesc.pkl'))
drugs_annot_moldesc.PubChemID = drugs_annot_moldesc.PubChemID.astype(str)

drugs_annot = pd.merge(drugs_annot_moldesc, drugs_annot_fam, on='PubChemID', how='left', sort=False)
assert (drugs_annot.PubChemID.astype(str) != drugs_tani.index.astype(str)).sum() == 0, 'order must be the same'
drugs_annot = drugs_annot.fillna("-")
logging.info(f'shape tanimoto matrix {drugs_tani.shape}, shape annotation {drugs_annot.shape}')
## save general in csv
drugs_tani.to_csv(os.path.join('../Results/data_per_dataset/full_data', f'tani_drugs.csv'), index = False, sep=";")
drugs_annot.to_csv(os.path.join('../Results/data_per_dataset/full_data', f'df_annot_drugs.csv'), index = False, sep=";")


# for Proteins
prot_rmsd = pd.read_pickle(os.path.join(PKLS_PATH, 'RMSD_full_matrix.pkl')) #
prot_annot = pd.read_pickle(os.path.join(PKLS_PATH, 'final_protein_annot.pkl'))
list_prots = list(set(prot_annot.UniprotID.tolist()).intersection(set(prot_rmsd.index.tolist())))
list_prots.sort()
prot_rmsd_final = prot_rmsd.loc[list_prots, list_prots]

prot_annot = prot_annot[prot_annot.UniprotID.isin(list_prots)].drop_duplicates()
prot_annot = prot_annot.set_index('UniprotID').loc[list_prots]
assert prot_rmsd_final.shape[0] ==  prot_annot.shape[0], 'shapes do not coincide'
assert (prot_rmsd_final.index != prot_annot.index).sum() == 0, 'index do not coincide'
logging.info(f'# prots {prot_rmsd_final.shape[0]}')
# save files in csv
prot_annot.to_csv(os.path.join('../Results/data_per_dataset/full_data', f'prot_annot.csv'), index = False, sep=";")
prot_rmsd_final.to_csv(os.path.join('../Results/data_per_dataset/full_data', f'prot_rmsd_final.csv'), index = False, sep=";")


# Load all DTIS in a dicionary
dict_dfs  = {
             'DrugBank': hf.get_dtis_drugbank(),
             'Davis_et_al': hf.get_dtis_davis(),
             'BindingDB': hf.get_bindingdb_dtis(),
             'BIOSNAP': hf.get_biosnap_dtis()
             }

dict_dfs.update(hf.get_dict_dtis_yamanishi())



### check if all folders exists if not create them
for key in list(dict_dfs.keys()):

    out_path = f'../Results/data_per_dataset/{key}'
    check_folder = os.path.isdir(out_path)
    if not check_folder:
        os.makedirs(out_path)
        logging.debug(f"created folder : {out_path}")
    else:
        logging.debug(f"{out_path} folder already exists.")


for key in list(dict_dfs.keys()):
    out_path = f'../Results/data_per_dataset/{key}'
    logging.info(f'Working in {key}')
    df = dict_dfs.get(key)
    drugs = df.PubChemID.astype(str).unique().tolist()
    prots = df.UniprotID.astype(str).unique().tolist()
    # save dataframe and list drugs & proteins
    df.to_csv(os.path.join(out_path, f'dtis_{key}.csv'), index = False, sep=";")
    np.savetxt(os.path.join(out_path, f'drugs_{key}.txt'), drugs, fmt='%s')
    np.savetxt(os.path.join(out_path, f'prots_{key}.txt'), prots, fmt='%s')
    # create tanimoto matrix
    df_drugs = drugs_tani.loc[drugs_tani.index.astype(str).isin(drugs), drugs_tani.columns.astype(str).isin(drugs)]
    df_drugs = df_drugs.iloc[~df_drugs.index.duplicated(keep='first'),~df_drugs.columns.duplicated(keep='first')]
    logging.debug(f' loosing: {set(drugs).difference(df_drugs.index.tolist())}')
    df_drugs.to_csv(os.path.join(out_path, f'drugs_sim_{key}.csv'), index = False, sep=";")
    assert (df_drugs.shape[0] <= len(drugs)) and (df_drugs.shape[1] <= len(drugs)) and (df_drugs.shape[0] == df_drugs.shape[1]), 'Not coinciding shapes in tanimoto matrix'
    logging.debug(df_drugs.head(3))
    # Annotations
    # do fill na with "-"
    df_drugs_annot = drugs_annot[drugs_annot.PubChemID.astype(str).isin(drugs)].drop_duplicates(subset='PubChemID')
    assert df_drugs_annot.shape[0] <= len(drugs), 'Not coinciding shapes in drug annotation data frame'
    assert (df_drugs_annot.PubChemID.astype(str) != df_drugs.index.astype(str)).sum() == 0, 'not coinciding order'
    df_drugs_annot = df_drugs_annot.fillna("-")
    df_drugs_annot.to_csv(os.path.join(out_path, f'drugs_annot_{key}.csv'), index = False, sep=";")
    logging.debug(df_drugs_annot.head(3))
    # check that index is the same !
    ####
    # rmsd protein matrix
    df_protein_rmsd = prot_rmsd.loc[prot_rmsd.index.isin(prots), prot_rmsd.columns.isin(prots)]
    assert df_protein_rmsd.shape[0] <=  len(prots), 'shape problem in rmsd matrix'
    assert df_protein_rmsd.index.shape ==  df_protein_rmsd.index.unique().shape, 'not unique index or column in rsmd martix'
    df_protein_rmsd.to_csv(os.path.join(out_path, f'prots_rmsd_{key}.csv'), index = False, sep=";")
    ####
    # protein annotations : NEED TO BE ORDERED!
    df_protein_annot = prot_annot[prot_annot.index.isin(prots)]
    assert df_protein_annot.shape[0] > 0, 'empty dataframe protein annotation data frame'
    # redo
    new_annot = pd.DataFrame(df_protein_rmsd.index, columns=['UniprotID'])
    new_annot = pd.merge(new_annot, df_protein_annot.reset_index(), on='UniprotID', how='left', sort=False)
    new_annot = new_annot.fillna("-")
    assert (new_annot.UniprotID.astype(str) != df_protein_rmsd.index.astype(str)).sum() == 0, 'not coinciding order'
    new_annot.to_csv(os.path.join(out_path, f'prots_annot_{key}.csv'), index = False, sep=";")
    print("")


