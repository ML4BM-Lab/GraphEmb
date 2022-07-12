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
import time
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import sys

# Using PubChemID Return all drugs in PubChemID
# Calculate similarity with Tanimoto
#  -> Fingerprint by default in RDKit
#  -> Family Annotations from ClassyFire

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)

##########


def get_SMILES_n_InCh_from_Pubchem_web_batch(drugs, size=500):
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
                inchikey =  [prop.get('value').get('sval') for prop in id.get('props') if prop.get('urn').get('label') == 'InChIKey'][0]
                if smile:
                    try:
                        mol1 = Chem.MolFromSmiles(str(smile))
                        fp1  = Chem.RDKFingerprint(mol1)
                    except:
                        logging.info(f'Error for pubchemid {cid}')
                        smile = None
                res.append((cid, smile, inchikey))
    return res


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




def get_drug_moldesc(smiles):
    mol = Chem.MolFromSmiles(smiles)
    log_p = Descriptors.MolLogP(mol)
    molwt = Descriptors.MolWt(mol)
    numhet = Descriptors.NumHeteroatoms(mol)
    return [log_p, molwt, numhet]





def get_drug_class(drugid, inkey, adapter):
    # search for classification
    retrieve_list = ['kingdom', 'superclass', 'class', 'subclass']
    #url = f'http://classyfire.wishartlab.com/entities/{inkey}.json'
    #response = requests.get(url)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)
    try:
        response = http.get(f'http://classyfire.wishartlab.com/entities/{inkey}.json')
    except:
        logging.info(f'connection did not work {sys.exc_info()[0]}')
    res = []
    if response.status_code == 200:
        json = response.json()
        clasif = [json.get(item, None).get('name', '-') for item in retrieve_list if json.get(item, None)]
        if len(clasif)<4:
            clasif.extend([None]*(4-len(clasif)))
        res.extend(clasif)        
    else:
        logging.info(f'status code {response.status_code} for {drugid} : {inkey}')
        results_.extend([None]*4)
    # Molecular descriptors
    res.insert(0, drugid)
    return res



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


df_smiles = pd.DataFrame(list_drug_smiles_inch, columns= ['drugid', 'smiles', 'inchkey']) 
df_smiles = df_smiles.drop_duplicates().dropna()

############################################ 
## TANIMOTO SCORE FOR CHEMICAL SIMILITUDE ##

# list_drugs = [str(code) for code, _,_ in list_drug_smiles_inch]
# list_smiles = [str(smiles) for _, smiles,_ in list_drug_smiles_inch] 
# list_inchi = [str(inchi) for _,_,inchi in list_drug_smiles_inch]
list_drugs = df_smiles.drugid.astype(str).tolist()
list_smiles = df_smiles.smiles.astype(str).tolist()
list_inchkey = df_smiles.inchkey.astype(str).tolist()

drug2smiles = dict(zip(df_smiles.drugid.astype(str), df_smiles.smiles.astype(str)))

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


df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns = list_drugs, index = list_drugs) 
logging.debug(f'Drug Similarity Matrix Shape {df_all_Sim_Tani.shape}')
df = df_all_Sim_Tani.astype(float)
# df.index = df.index.astype(str)
# df.columns = df.columns.astype(str)
df.to_pickle('../Data/pkls/drugs_tani_full.pkl')
#df = pd.read_pickle('../Data/pkls/drugs_tani_full.pkl')


###########################
####### Annotations #######

## create dict smiles

# annotations_moldesc = []
# for i in tqdm(range(len(list_smiles))):
#     line = get_drug_moldesc(list_smiles[i])
#     line.insert(0, list_drugs[i])
#     annotations_moldesc.append(line)

# try instead
annotations_moldesc = []
for drug in tqdm(df.index):    
    sml = drug2smiles.get(drug, None)
    line = [drug]
    line.extend(get_drug_moldesc(sml))
    annotations_moldesc.append(line)

annotations_moldesc
df_annot_moldesc = pd.DataFrame.from_records(annotations_moldesc, columns=['PubChemID',  'log_p', 'molwt', 'numhet'])
df_annot_moldesc.to_pickle('../Data/pkls/annot_moldesc.pkl')


## 
retry_strategy = Retry(
    total=3,
    status_forcelist = [429], # , 500, 502, 503, 504
    method_whitelist=["HEAD", "GET", "OPTIONS"]
)
adapter = HTTPAdapter(max_retries=retry_strategy)

# create dict inchkey
annotations_family = []
for i in tqdm(range(len(list_drugs)), desc='Annotating drugs'):
    drugid = list_drugs[i]
    smiles = list_smiles[i]
    inkey = list_inchkey[i]
    results_ =  get_drug_class(drugid, inkey, adapter)
    annotations_family.append(tuple(results_))


##
df_annot_cols = ['PubChemID', 'kingdom', 'superclass', 'class', 'subclass']
df_annot_fam = pd.DataFrame.from_records(annotations_family, columns=df_annot_cols)
df_annot_fam.to_pickle('../Data/pkls/df_annot_fam.pkl')



Annot_Drugs = pd.DataFrame(df.index.tolist(), columns=['PubChemID'])

final_Annot_Drugs = pd.merge(Annot_Drugs, df_annot.drop_duplicates(), on='PubChemID', how='left').fillna('-')
assert (df.index != final_Annot_Drugs.PubChemID).sum() == 0, 'order do not coincide'
logging.info(f'shape tanimoto matrix {df.shape}, shape annotation {final_Annot_Drugs.shape}')

final_Annot_Drugs.to_pickle('../Data/pkls/drugs_all_annotation.pkl')

#df_annot = pd.read_pickle('../Data/pkls/drugs_annot_full.pkl')

## MOL DESCRIPTORS

# Number Heteroatoms

# mol_descrt_list = []
# for drug, smiles in zip(list_drugs, list_smiles):
#     m = Chem.MolFromSmiles(smiles)
#     log_p = Descriptors.MolLogP(m)
#     molwt = Descriptors.MolWt(m)
#     numhet = Descriptors.NumHeteroatoms(m)
#     res = (drug, log_p, molwt, numhet)
#     mol_descrt_list.append(res)


# cols_mol_descript = ['PubChemID', 'log_p', 'molwt', 'numhet']
# df_mol_descr = pd.DataFrame.from_records(mol_descrt_list, columns=cols_mol_descript)

# # join 2 dataframes

# all_anot_drugs = df_annot.merge(df_mol_descr, on='PubChemID')

# # frames = [df_annot, df_mol_descr]
# # all_anot_drugs = pd.concat(frames, axis =1)
# # test_shape = (df_annot.shape[0], df_annot.shape[1]+df_mol_descr.shape[1])

# assert all_anot_drugs.shape[0] == df.shape[0], 'Shapes do not coincide'

# all_anot_drugs.to_pickle('../Data/pkls/drugs_all_annotation.pkl')

# logging.debug('0')
