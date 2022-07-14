
import os
import pandas as pd
import numpy as np
import helper_functions as hf
import logging
import requests
from tqdm import tqdm
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import sys
import time
from rdkit import Chem
import re
from rdkit.Chem import Descriptors

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)



####

def get_key(ecode, kws):
    global class_list
    class_list = np.loadtxt('../Data/keys.txt', dtype='str', delimiter='\n').tolist()
    class_list = [key.lower() for key in class_list]
    logging.debug("""keywords from: https://www.uniprot.org/keywords/KW-9992; 
                    deleting enzyme keywords because using EC Code for those""")
    if ecode:
        mol_funct = 'enzyme'
    elif any(key in class_list for key in kws):
        mol_funct = [word for word in class_list if word in kws][0]
    else:
        mol_funct = 'notAnnotated' 
    return mol_funct


def get_json_info_uni(data):
    check_uni_id = data['primaryAccession']
    length = int(data['sequence']['length'])
    ecode = None
    if data.get('proteinDescription').get('recommendedName'):
        recname = data.get('proteinDescription').get('recommendedName')
        if recname.get('ecNumbers'):
            ecode = recname.get('ecNumbers')[0].get('value')
    mol_funct = None
    if data.get('keywords'):
        dic_molfunc = [line for line in data.get('keywords') if line.get('category') == 'Molecular function']
        keywords = [line.get('name').lower() for line in dic_molfunc]
        mol_funct = get_key(ecode, keywords) 
    res = (check_uni_id, length, ecode, mol_funct)
    return res


def get_prot_function(proteins):
    #list_proteins =  list_proteins_struct # np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()
    nprots = len(proteins)
    if nprots < 500:
        n = int(len(proteins)/2)
    else:
        n = 500
    #
    list_proteins = [proteins[i:i + n] for i in range(0, len(proteins), n)]
    #
    full_info= []
    for lts in tqdm(list_proteins, desc='Retrieving uniprot annotation'):
        unilist = ','.join(lts)
        r = requests.get(f'https://rest.uniprot.org/uniprotkb/accessions?accessions={unilist}')
        jsons = r.json()['results']
        for data in jsons:
            res = get_json_info_uni(data)
            full_info.append(res)
    #
    # requests lost conection using time just to be ok
    time.sleep(2)
    #
    if len([(id, func) for id,_,_,func in full_info]) == nprots:
        logging.info('all proteis as principal accesion id')
        protein2function = dict([(id, func) for id,_,_,func in full_info])
    else:
        logging.info('checking in secondary accesion numbers...')
        retrn_prots = [id for id,_,_,_ in full_info]
        lost_proteins = list(set(proteins) - set(retrn_prots))
        logging.info(f'{len(lost_proteins)} proteins changed un uniprot, accesing secondary')
        # those can be queried with the secondary accesion number 
        for prot in tqdm(lost_proteins, desc='Secondary Acc'):
            r = requests.get(f'https://rest.uniprot.org/uniprotkb/search?&query=sec_acc:{prot}')
            data = r.json()['results'][0]
            res = get_json_info_uni(data)
            full_info.append(res)
        retrn_prots = [id for id,_,_,_ in full_info]
        lost_proteins = list(set(proteins) - set(retrn_prots))
        logging.info(f'loosing {len(lost_proteins)} proteins')
        protein2function = dict([(id, func) for id,_,_,func in full_info])
    return protein2function

###

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
        res.extend([None]*4)
    # Molecular descriptors
    res.insert(0, drugid)
    return res


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


def get_drug_superclass(drugs):
    drug2inchkey = dict([(str(drugid), inchkey) for drugid,_,inchkey in get_SMILES_n_InCh_from_Pubchem_web_batch(drugs)])
    ## 
    retry_strategy = Retry(
        total=3,
        status_forcelist = [429], # , 500, 502, 503, 504
        method_whitelist=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    # create dict inchkey
    annotations_family = []
    for drug in tqdm(drugs, desc='Annotating drugs'):
        inkey = drug2inchkey.get(drug, None)
        results_ =  get_drug_class(drug, inkey, adapter)
        annotations_family.append(tuple(results_))
    #
    drug2superclass = dict([(id, superclass) for id,_,superclass,_,_ in annotations_family ])
    return drug2superclass


def get_drug_subclass(drugs):
    drug2inchkey = dict([(str(drugid), inchkey) for drugid,_,inchkey in get_SMILES_n_InCh_from_Pubchem_web_batch(drugs)])
    ## 
    retry_strategy = Retry(
        total=3,
        status_forcelist = [429], # , 500, 502, 503, 504
        method_whitelist=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    # create dict inchkey
    annotations_family = []
    for drug in tqdm(drugs, desc='Annotating drugs'):
        inkey = drug2inchkey.get(drug, None)
        results_ =  get_drug_class(drug, inkey, adapter)
        annotations_family.append(tuple(results_))
    #
    drug2subclass = dict([(id, subclass) for id,_,_,_,subclass in annotations_family ])
    return drug2subclass

############################################################
###### START ###

logging.info('Loading dtis')
# Load all DTIS in a dicionary
dict_dfs  = {
             'DrugBank': hf.get_dtis_drugbank(),
             'Davis_et_al': hf.get_dtis_davis(),
             'BindingDB': hf.get_bindingdb_dtis(),
             'BIOSNAP': hf.get_biosnap_dtis()
             }

dict_dfs.update(hf.get_dict_dtis_yamanishi())


### specific dataset
# load
for key in dict_dfs.keys():
    out_path = f'../Results/data_per_dataset/{key}'
    logging.info(f'Working in {key}')
    df = dict_dfs.get(key)
    # get dicts
    # before only retrieved for crystalized now for each dataset for all
    proteins = df.UniprotID.tolist()
    protein2function = get_prot_function(proteins)
    #
    #### drug annotation
    drugs = df.PubChemID.tolist()
    drug2superclass = get_drug_superclass(drugs)
    # drug2subclass = get_drug_subclass(drugs)
    #
    # append to df
    df['prot_func'] = df.UniprotID.map(protein2function)
    df['drug_fam'] = df.PubChemID.map(drug2superclass)
    df_grouped = df.groupby(['drug_fam', 'prot_func'])['drug_fam'].count()
    data_plot = pd.DataFrame([(a[0], a[1],c) for a,c in zip(df_grouped.index,df_grouped.tolist())], columns = ['source', 'target', 'count'])
    data_plot.to_pickle(f'../Results/sankeys/data_grouped_{key}.pkl')
    data_plot.to_csv(f'../Results/sankeys/data_grouped_{key}.csv', sep=';', index=False)
    logging.info('--')

