import os
import pandas as pd
import logging
from tqdm import tqdm
import helper_functions as hf
import requests
from tqdm import tqdm
import numpy as np
import requests
from itertools import repeat
from rdkit import Chem
import re
from rdkit.Chem import Descriptors

# Using PubChemID Return all drugs in PubChemID
# Calculate similarity with Tanimoto
#  -> Fingerprint by default in RDKit
#  -> Family Annotations from ClassyFire

logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)

##########

def get_SMILES_n_InCh_from_Pubchem_web_batch(drugs):
    smiles = []
    # split the list in chunks of 100 to make requests
    size = 100
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
        elif response.status_code == 404:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/sid/{chunk}/json"
            response = requests.get(url)
            jsons = response.json()
        for id in jsons.get("PC_Compounds"):
            cid = id.get('id').get('id').get('cid')
            smile = [prop.get('value').get('sval') for prop in id.get('props') if prop.get('urn').get('label') == 'SMILES' and prop.get('urn').get('name') == 'Canonical']
            inchikey =  [prop.get('value').get('sval') for prop in id.get('props') if prop.get('urn').get('label') == 'InChIKey']
            if smile:
                try:
                    mol1 = Chem.MolFromSmiles(str(smile[0]))
                    fp1  = Chem.RDKFingerprint(mol1)
                    smiles.append((cid, smile[0], inchikey[0]))
                except:
                    logging.info(f'Error for pubchemid {cid}')
    return smiles


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


#################################################### 
################ Load Drugs & SMILES ###############
drugs_drugbank = get_drugs_drugbank() # only drugbank
drugs_davis = get_drugs_Davis()
drugs_binding = get_drugs_binding()
drugs_biosnap = get_drugs_biosnap()
drugs_yamanishi = get_drugs_yamanishi()

list_drugs_ = drugs_drugbank + drugs_davis + drugs_binding + drugs_biosnap + drugs_yamanishi
list_drugs_ = [str(drug) for drug in list_drugs_]
all_drugs = list(set(list_drugs_))
#all_drugs_2 = list(set([str(drug) for drug in all_drugs]))

logging.info(f'Total number of drugs: {len(all_drugs)}')
list_drug_smiles_inch = get_SMILES_n_InCh_from_Pubchem_web_batch(all_drugs) 


############################################ 
## TANIMOTO SCORE FOR CHEMICAL SIMILITUDE ##

list_drugs = [code for code, _, _ in list_drug_smiles_inch]
list_smiles = [smiles for _, smiles,_ in list_drug_smiles_inch] 
list_inchi = [inchi for _,_,inchi in list_drug_smiles_inch]

# Loop for creating Similarity Matrix
all_Sim_Tani = []
dic = {}
for drug in tqdm(range(len(list_drugs)), desc='Retrieving Pairwise Tanimoto for drugs', position=0, leave=True): 
    id_drug_1 = list_drugs[drug]
    smiles_drug_1 = list_smiles[list_drugs.index(id_drug_1)]
    sim_for_drug = []
    tmp = []
    tmp.extend(repeat(smiles_drug_1, len(list_drugs))) 
    for j in range(len(tmp)):
        result = hf.get_pairwise_tanimoto(tmp[j], list_smiles[j], dic)
        sim_for_drug.append(result)
    all_Sim_Tani.append(sim_for_drug)


df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= list_drugs, index = list_drugs) 
logging.debug(f'Drug Similarity Matrix Shape {df_all_Sim_Tani.shape}')
df = df_all_Sim_Tani.astype(float)
df.index = df.index.astype(str)
df.columns = df.columns.astype(str)
df.to_pickle('../Data/pkls/drugs_tani_full.pkl')

###########################
####### Annotations #######

### CLASS (CLASSYFIRE)
annotations = []
retrieve = ['kingdom', 'superclass', 'class', 'subclass']
for i in tqdm(range(len(list_inchi)), desc='Annotating drugs'):
    elmt = list_inchi[i]
    url = f'http://classyfire.wishartlab.com/entities/{elmt}.json'
    response = requests.get(url)
    results_ = []
    if response.status_code == 200:
        json = response.json()
        for item in retrieve:
            t = json.get(item, None)
            if t: 
                result = t.get('name', None)
                results_.append(result)
    results_.insert(0, list_drugs[i])
    annotations.append(tuple(results_))


df_annot = pd.DataFrame.from_records(annotations, columns=['PubChemID', 'kingdom', 'superclass', 'class', 'subclass'])
df_annot.to_pickle('../Data/pkls/drugs_annot_full.pkl')


## MOL DESCRIPTORS

# Number Heteroatoms

smile_test = list_smiles[2]

mol_descrt_list = []
for drug, smiles in zip(list_drugs, list_smiles):
    m = Chem.MolFromSmiles(smiles)
    log_p = Descriptors.MolLogP(m)
    molwt = Descriptors.MolWt(m)
    numhet = Descriptors.NumHeteroatoms(m)
    res = (drug, log_p, molwt, numhet)
    mol_descrt_list.append(res)


cols_mol_descript = ['PubChemID', 'log_p', 'molwt', 'numhet']
df_mol_descr = pd.DataFrame.from_records(mol_descrt_list, columns=cols_mol_descript)

# join 2 dataframes
frames = [df_annot, df_mol_descr]
all_anot_drugs = pd.concat(frames, axis =1)
test_shape = (df_annot.shape[0], df_annot.shape[1]+df_mol_descr.shape[1])
assert all_anot_drugs.shape == test_shape, 'Shapes do not coincide'

all_anot_drugs.to_pickle('../Data/pkls/drugs_all_annotation.pkl')

logging.debug('0')
