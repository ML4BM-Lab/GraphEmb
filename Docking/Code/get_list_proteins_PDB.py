import os
import pandas as pd
import logging
from tqdm import tqdm
import helper_functions as hf
import requests
from tqdm import tqdm
#import glob 
#from collections import Counter
import numpy as np
#from Bio.PDB import PDBParser
#from Bio.PDB.MMCIFParser import MMCIFParser
#import gzip


### functions
def get_unis_yamanishi(all_prot):
    GEN_PATH = '../../DB/Data/'
    PATH_E = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt')
    PATH_NR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt')
    PATH_GPCR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt')
    PATH_IC = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt')
    list_yam = [PATH_E, PATH_NR, PATH_GPCR, PATH_IC]
    for yam_db in list_yam:
        print(yam_db)
        # yam
        colnames_data = ['Kegg_ID', 'Gene']
        df = pd.read_csv(yam_db, header = None, names = colnames_data, index_col=False, sep='\t')
        df['Gene'] = df['Gene'].map(lambda x: x[:3] + ":" + x[3:])
        # change from hsa: to Uniprot ---> wget 
        logging.info('Loading hsa2uni dictionary...')
        hsa2uni = hf.get_dict_hsa2uni()
        # repeated entries, but diff uniprots have the same sequence (so it's safe)
        df['UniprotID'] = df['Gene'].map(hsa2uni)
        logging.debug(f'{(df.head(2))}')
        df = df.dropna()
        prots = df.UniprotID.unique().tolist()
        all_prot.extend(prots)
        all_prot = list(set(all_prot))
        logging.debug(len(all_prot))
    return all_prot

#  
def get_unis_drugbank(all_prot):
    GEN_PATH = '../../DB/Data/'
    PATH_DRUGBANK = os.path.join(GEN_PATH, 'DrugBank/DrugBank_DTIs.tsv')
    df = pd.read_csv(PATH_DRUGBANK, sep='\t') 
    df.columns = ['Drug', 'Protein']
    df = df.drop_duplicates()
    prots = df.Protein.unique().tolist()
    all_prot.extend(prots)
    all_prot = list(set(all_prot))
    logging.debug(len(all_prot))
    return all_prot

# 
def get_unis_biosnap(all_prot):
    GEN_PATH = '../../DB/Data/'
    PATH_BIOSNAP = os.path.join(GEN_PATH, 'BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv')
    df = pd.read_csv(PATH_BIOSNAP, sep='\t', comment='#', header=None) 
    df.columns = ['Drug', 'Protein']
    df = df.drop_duplicates()
    prots = df.Protein.unique().tolist()
    #len(all_prot)
    all_prot.extend(prots)
    all_prot = list(set(all_prot))
    logging.debug(len(all_prot))
    return all_prot

# 
def get_unis_binding(all_prot):
    GEN_PATH = '../../DB/Data/'
    PATH_BINDING =  os.path.join(GEN_PATH,'BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv')
    #
    bindingdb = pd.read_csv(PATH_BINDING, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Y'])
    bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID', 'Target Sequence': 'Sequence'}, axis=1)
    # threshold
    mod_bind = bindingdb.copy()
    threshold = 30
    mod_bind['Label'] = [1 if x < threshold else 0 for x in mod_bind['Y']]
    mod_bind = mod_bind.drop(columns='Y')
    logging.debug(f'Shape before removing Label ==0 {mod_bind.shape}')
    mod_bind = mod_bind.drop(mod_bind[mod_bind.Label == 0].index)
    logging.debug(f'Shape after removing 0s {mod_bind.shape}')
    mod_bind = mod_bind.drop(columns='Label')
    mod_bind = mod_bind.dropna()
    prots = mod_bind.UniprotID.unique().tolist()
    all_prot.extend(prots)
    all_prot = list(set(all_prot))
    logging.debug(len(all_prot))
    return all_prot

# 
def get_unis_davis(all_prot):
    GEN_PATH = '../../DB/Data/'
    PATH_DAVIS = os.path.join(GEN_PATH, 'Davis_et_al/tdc_package_preprocessing/DAVIS_et_al_w_labels.tsv')
    davis = pd.read_csv(PATH_DAVIS, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Label'])
    davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
    mod_davis = davis.copy()
    # then drop those with 0
    mod_davis = mod_davis.drop(mod_davis[mod_davis.Label == 0].index)
    mod_davis = mod_davis.drop(columns=['PubChemID','SMILES','Label'])
    ## PROTEIN IDENTIFIERS
    data_path_biomart = os.path.join(GEN_PATH, 'cross_side_information_DB/bioMART/mart_export_expanded.txt')
    biomart = pd.read_csv(data_path_biomart, sep='\t', usecols=['UniProtKB Gene Name symbol', 'UniProtKB Gene Name ID'])
    biomart = biomart.dropna().drop_duplicates()
    genename2geneid = dict(zip(biomart['UniProtKB Gene Name symbol'].tolist(), biomart['UniProtKB Gene Name ID'].tolist()))
    logging.debug('mapping...')
    mod_davis['UniprotID'] = mod_davis['GeneName'].map(genename2geneid)
    mod_davis = mod_davis.dropna().drop_duplicates()
    prots = mod_davis.UniprotID.unique().tolist()
    all_prot.extend(prots)
    all_prot = list(set(all_prot))
    logging.debug(len(all_prot))
    return all_prot
#
## all unis join
def get_all_unis():
    # get list of all uniprots
    all_prot = []
    all_prot = get_unis_yamanishi(all_prot)
    all_prot = get_unis_drugbank(all_prot)
    all_prot = get_unis_biosnap(all_prot)
    all_prot = get_unis_binding(all_prot)
    all_prot = get_unis_davis(all_prot)
    all_prot = list(set(all_prot)) # safety check
    return all_prot


def get_df_info_uni2pdb(all_prot):
    df_uni2pdb = pd.DataFrame(columns=['UniprotID', 'PDB', 'Method', 'Resolution', 'chain-resid'])
    for uniprotid in tqdm(all_prot):
        r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
        a = r.text.split('\n')
        check_pdb = any('PDB;' in string for string in a)
        if check_pdb:
            pdb_info = [line for line in a if ' PDB; ' in line]
            for s in pdb_info:
                info = s.split(';')
                info = [nfo.strip() for nfo in info]
                df_uni2pdb.loc[len(df_uni2pdb.index)] = [uniprotid, info[1], info[2], info[3], info[4]]
    return df_uni2pdb



######################
logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)
######################

all_prot = get_all_unis()
logging.info(f'Number of all proteins: {len(all_prot)}') # 6167

df_uni2pdb = get_df_info_uni2pdb(all_prot)
df_uni2pdb = df_uni2pdb.drop(df_uni2pdb[df_uni2pdb.Method != 'X-ray'].index)
df_uni2pdb = df_uni2pdb.drop(df_uni2pdb[df_uni2pdb.Resolution == '-'].index)
df_uni2pdb.Resolution = df_uni2pdb.Resolution.str.strip(' A').astype(float)
df_uni2pdb = df_uni2pdb.drop(columns='Method') # only X-Ray

logging.debug(df_uni2pdb.head(6))

df_uni2pdb.sort_values(by='Resolution', inplace=True, ignore_index=True) # sorting values by resolution, better will be on top of the dataframe
df_uni2pdb.drop_duplicates(subset='UniprotID', inplace=True, ignore_index=True) # remove uniprots duplicated, will keep only the first appearing (sorted before)
df_uni2pdb.sort_values(by=['UniprotID'], inplace=True, ignore_index=True) # Now sorting again by uniprotID alphabet
res_uni2pdb =  df_uni2pdb.drop(df_uni2pdb[df_uni2pdb.Resolution > 2].index)

res_uni2pdb.to_pickle('pdb_data/res_uni2pdb.pkl')
# res_uni2pdb = pd.read_pickle('../Data/res_uni2pdb.pkl')

res_uni2pdb.shape
logging.info(f'Number of proteins (uniprot ID) with available PDB with res<2 A: {len(res_uni2pdb.UniprotID.unique())}')

# save and download again!!!
unique_pdbs =  res_uni2pdb.PDB.unique().tolist() #batch download
file_list = 'list_download_pdbs.txt'
np.savetxt(file_list, unique_pdbs , newline='', fmt='%s,')
