import os
import logging
import numpy as np
import requests
import re
import pandas as pd
from tqdm import tqdm
import helper_functions as hf
from rdkit.Chem import Descriptors
from rdkit import Chem
import seaborn as sns
import matplotlib.pyplot as plt

def get_cid_from_sid(drugs):
    '''
    returns a dict
    '''
    cids = []
    size = 100
    drugs = [str(drug) for drug in drugs]
    split_list = lambda big_list, x: [
        big_list[i : i + x] for i in range(0, len(big_list), x)
    ]
    drug_chunks = split_list(drugs, size)
    for chunk in tqdm(drug_chunks, desc='Requesting SIDs in PubChem'):
        chunk = ",".join(chunk)
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{chunk}/json'
        response = requests.get(url)
        if response.status_code == 200:
            jsons = response.json()
        for id in jsons.get("PC_Substances"):
            try:
                id = id.get('compound', None)
                cid = id[1].get('id').get('id').get('cid')
            except:
                print(f'Differnt format for drug')
                cid = None
            cids.append(cid)
    return dict(zip(drugs, cids))


def get_Kegg_pubchem_dict():
    PUBCHEM_PATTERN = re.compile(r"(?<=pubchem:)[\d]+")
    KEGG_PATTERN = re.compile(r"(?<=dr:)[\w\d]+")
    r = requests.get(f"http://rest.kegg.jp/conv/drug/pubchem")
    if r.status_code == 200:
        dictionary = r.text.strip().split("\n")
        dictionary = [
            (KEGG_PATTERN.search(entry).group(), PUBCHEM_PATTERN.search(entry).group())
            for entry in dictionary
        ]
        return dict(dictionary)

def get_drugs_yamanishi():
    # defining paths
    GEN_PATH = '../../DB/Data/'
    PATH_E = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt')
    PATH_NR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt')
    PATH_GPCR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt')
    PATH_IC = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt')
    list_yam = [PATH_E, PATH_NR, PATH_GPCR, PATH_IC]
    # Retrieving all drugs
    drugs_yam_kegg = []
    for yam_db in list_yam:
        yam_db
        colnames_data = ['Kegg_ID', 'Gene']
        df = pd.read_csv(yam_db, header = None, names = colnames_data, index_col=False, sep='\t')
        drugs_yam_kegg.extend(df.Kegg_ID.unique().tolist())
    drugs_yam_kegg = list(set(drugs_yam_kegg))
    logging.info(f'# of drugs with Kegg ID {len(drugs_yam_kegg)}')
    kegg2sid = get_Kegg_pubchem_dict()
    drugs_yam_sid = [kegg2sid.get(drug,None) for drug in drugs_yam_kegg]
    sid2cid = get_cid_from_sid(drugs_yam_sid)
    drugs_yam_pc = [sid2cid.get(drug, None) for drug in drugs_yam_sid]
    return drugs_yam_pc

def get_drugs_yamanishi_sep():
    # defining paths
    GEN_PATH = '../../DB/Data/'
    PATH_E = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt')
    PATH_NR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt')
    PATH_GPCR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt')
    PATH_IC = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt')

    list_yam = [PATH_E, PATH_NR, PATH_GPCR, PATH_IC]
    dict_yamanishis = dict(zip(['e', 'nr', 'gpcr', 'ic'], list_yam))

    dict_yamanishis_drugs_pc = {}

    #for yam_db in list_yam:
    for name_db in list(dict_yamanishis.keys()):
        yam_db = dict_yamanishis[name_db]

        colnames_data = ['Kegg_ID', 'Gene']
        df = pd.read_csv(yam_db, header = None, names = colnames_data, index_col=False, sep='\t')
        drugs_kegg = df.Kegg_ID.unique().tolist()

        logging.info(f'For {name_db} # of drugs with Kegg ID {len(drugs_kegg)}')
        kegg2sid = get_Kegg_pubchem_dict()

        drugs_yam_sid = [kegg2sid.get(drug,None) for drug in drugs_kegg]

        sid2cid = get_cid_from_sid(drugs_yam_sid)

        drugs_pc = [sid2cid.get(drug, None) for drug in drugs_yam_sid]

        # filter here nones:
        drugs_pc = [str(drug) for drug in drugs_pc if drug is not None] # check

        dict_yamanishis_drugs_pc[name_db] = drugs_pc

    return dict_yamanishis_drugs_pc


# Bindingdb 
def get_drugs_binding():
    db_file_path = '../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
    bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Y'])
    bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID', 'Target Sequence': 'Sequence'}, axis=1)
    mod_bind = bindingdb.copy()
    threshold = 30
    mod_bind['Label'] = [1 if x < threshold else 0 for x in mod_bind['Y']]
    mod_bind = mod_bind.drop(columns='Y')
    mod_bind = mod_bind.drop(mod_bind[mod_bind.Label == 0].index)
    mod_bind = mod_bind.drop(columns='Label')
    mod_bind.loc[:, 'PubChemID'] = mod_bind.loc[:, 'PubChemID'].astype(int).astype(str) # not float
    drugs_binding = mod_bind.PubChemID.unique().tolist()
    return drugs_binding

## biosnap
def get_drugs_biosnap():
    dti_file_path = '../../DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
    dti = pd.read_csv(os.path.join(os.getcwd(), dti_file_path), sep='\t', comment='#', header=None) 
    dti.columns = ['Drug', 'Protein']
    drugs = dti.Drug.unique().tolist()
    # get dict and inverse it
    dic_cid_dbid =  hf.pubchem_to_drugbankid()
    dic_dbid_cid = {v: k for k, v in dic_cid_dbid.items()}
    drugs_pc = [dic_dbid_cid.get(i, None) for i in drugs] # 
    drugs_pc = list(filter(None, drugs_pc)) 
    return drugs_pc

# Davis
def get_drugs_Davis():
    db_file_path = '../../DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al_w_labels.tsv'
    davis = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Label'])
    davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
    logging.debug(davis.head(2))
    mod_davis = davis.copy()
    # we only need the positive edges for writing a edge list
    # then drop those with 0
    logging.debug('Drop null int for edge list')
    logging.debug(f'Shape before removing Label ==0 {mod_davis.shape}')
    mod_davis = mod_davis.drop(mod_davis[mod_davis.Label == 0].index)
    logging.debug(f'Shape after removing 0s {mod_davis.shape}')
    mod_davis = mod_davis.drop(columns='Label')
    drugs_davis = mod_davis.PubChemID.unique().tolist()
    return drugs_davis

### Drugbank DB
def get_drugs_drugbank():
    PATH_DB = '../../DB/Data/DrugBank/DrugBank_DTIs.tsv'
    dti = pd.read_csv(os.path.join(os.getcwd(), PATH_DB), sep='\t') 
    dti.columns = ['Drug', 'Protein']
    drugs = dti.Drug.unique().tolist() # in DrugBankID
    logging.debug(f'drugbank drugs {len(drugs)}')
    # to PubChemID
    dic_cid_dbid =  hf.pubchem_to_drugbankid()
    dic_dbid_cid = {v: k for k, v in dic_cid_dbid.items()}
    drugs_pc = [dic_dbid_cid.get(i, None) for i in drugs] # 
    drugs_pc = list(filter(None, drugs_pc)) 
    return drugs_pc


def get_SMILES_Pubchem_web_batch(drugs, size=500):
    res = []
    # split the list in chunks of 100 to make requests
    drugs = [str(drug) for drug in drugs]
    split_list = lambda big_list, x: [
        big_list[i : i + x] for i in range(0, len(big_list), x)
    ]
    drug_chunks = split_list(drugs, size)
    for chunk in tqdm(drug_chunks, desc='Requesting SMILES to PubChem'):
        chunk = ",".join(chunk)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{chunk}/json"
        response = requests.get(url)
        if response.status_code == 200:
            jsons = response.json()
            for id in jsons.get("PC_Compounds"):
                cid, smile, inchikey = None, None, None
                cid = id.get('id').get('id').get('cid')
                smile = [prop.get('value').get('sval') for prop in id.get('props') if prop.get('urn').get('label') == 'SMILES' and prop.get('urn').get('name') == 'Canonical'][0]
                #inchikey =  [prop.get('value').get('sval') for prop in id.get('props') if prop.get('urn').get('label') == 'InChIKey'][0]
                if smile:
                    try:
                        mol1 = Chem.MolFromSmiles(str(smile))
                        fp1  = Chem.RDKFingerprint(mol1)
                        res.append((cid, smile))
                    except:
                        logging.info(f'Error for pubchemid {cid}')
                        smile = None
                        print('error')    
        else:
            print(response)                
    return res




# drugs_pubchem = dict_datasets.get('gpcr')

# list_drug_smiles_inch = get_SMILES_Pubchem_web_batch(drugs_pubchem, size=10) 

# drugs = [str(drug) for drug in drugs_pubchem]
# split_list = lambda big_list, x: [
#     big_list[i : i + x] for i in range(0, len(big_list), x)
# ]
# drug_chunks = split_list(drugs, 10)

# a = list(filter(lambda item: item is not 'None', drugs_pubchem))
# len(drugs_pubchem)
#################################################### 
################ Load Drugs & SMILES ###############

# there are 'None' n Yamanishi!! 
def calculate_fp(list_datasets):

    df = pd.DataFrame()

    for database in list_datasets:

        print(database)

        drugs_pubchem = dict_datasets.get(database)

        drugs_pubchem = list(filter(lambda item: item is not 'None', drugs_pubchem))

        list_drug_smiles_inch = get_SMILES_Pubchem_web_batch(drugs_pubchem, 10) 

        drugs_smiles = [j for _,j in list_drug_smiles_inch]
        drug_mol = [Chem.MolFromSmiles(smiles) for smiles in drugs_smiles]

        df[database] = pd.Series(drug_mol)
    
    return df


def df_moldesc(df_fp, mol_descr):
    
    df = df_fp.copy()

    dict_moldesc = {'mw': Descriptors.MolWt,
                    'logp':  Descriptors.MolLogP,
                    'HeavyAtomCount': Descriptors.HeavyAtomCount,
                    'NHOHCount': Descriptors.NHOHCount,
                    'NOCount': Descriptors.NOCount,
                    'NumAliphaticRings': Descriptors.NumAliphaticRings,
                    'NumAromaticRings': Descriptors.NumAromaticRings,
                    'NumHAcceptors': Descriptors.NumHAcceptors,
                    'NumHDonors': Descriptors.NumHDonors,
                    'NumHeteroatoms': Descriptors.NumHeteroatoms,
                    'NumRotatableBonds': Descriptors.NumRotatableBonds,
                    'RingCount': Descriptors.RingCount,
                    'NumValenceElectrons': Descriptors.NumValenceElectrons,
                    #'NumAmideBonds': Descriptors.NumAmideBonds,
                    'TPSA': Descriptors.TPSA
                    }
    
    f = dict_moldesc.get(mol_descr)

    for dataset in list_datasets: 
        df[dataset] = df[dataset].apply(lambda a: f(a) if pd.notnull(a) else a)
    
    return df



drugs_drugbank = get_drugs_drugbank() 
drugs_davis = get_drugs_Davis()
drugs_binding = get_drugs_binding()
drugs_biosnap = get_drugs_biosnap()
dic_drugs_yamanishi = get_drugs_yamanishi_sep() 


dict_datasets = {'drugbank': drugs_drugbank,
                 'biosnap': drugs_biosnap,
                 'bindingdb': drugs_binding,
                 'davis': drugs_davis
                 }

dict_datasets.update(dic_drugs_yamanishi)


list_datasets = list(dict_datasets.keys())

df_fp = calculate_fp(list_datasets)


mol_descr = 'mw'

list_md = [ 'mw', 'logp', 'HeavyAtomCount', 'NHOHCount', 'NOCount', 
           'NumAliphaticRings', 'NumAromaticRings', 'NumHAcceptors', 
           'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'RingCount',
           'NumValenceElectrons', 'NumAmideBonds', 'TPSA'
           ]


for mol_descr in list_md:

    df = df_moldesc(df_fp, mol_descr)

    plt.clf()
    sns.set_palette("Set2")
    plt.figure(figsize=(8, 8))
    plt.title(mol_descr)
    #plt.ylim(0, 2_000)
    #plt.xlim(0, 2_000)
    sns.boxplot(data=df)
    #sns.histplot(data=df, kde=True, common_norm=False, alpha=0.25) #displot
    plt.savefig(f'../Results/moldesc_boxplots/boxplot_{mol_descr}.pdf', dpi=330 ,bbox_inches='tight',  pad_inches = 0.25)

