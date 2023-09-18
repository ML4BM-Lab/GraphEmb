import os, sys, uuid
import argparse
import logging
import numpy as np
import pandas as pd
import json
import requests
from tqdm import tqdm
from re import search
import time
import xml.etree.ElementTree as ET
import pubchempy as pcp
import subprocess as sp
from pubchempy import Compound
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
import random
import re



'''
from importlib import reload
reload(hf)
'''

# Databases Paths
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

# funct
def get_cid(compound_name):
    '''
    Not used (?)
    This function returns the compound id (PubChem) given a name
    '''
    try:
        comp_id = pcp.get_compounds(str(compound_name), 'name')[0].cid 
        time.sleep(2) # increase if error
    except IndexError:
        comp_id = np.nan
    return comp_id

def get_compound_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

def pubchem_to_drugbankid():
    '''
    Parse  DrugBank DB xml and 
    create a dictionary from PubChemID to DrugBankID
    '''
    logging.debug('Reading DrugBank xml database')
    tree = ET.parse('../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    logging.debug('Succesfully read!')
    root = tree.getroot()
    dbids = []
    pubchemids = []
    for drug_entry in tqdm(root, desc='Retrieving Drugbank  & PubChem ID', position=0, leave=True):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        pubchem_id = np.nan # to not repeat id if it does not appear later
        for props in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
            for prop in props:
                if(prop.text == 'PubChem Compound'): 
                    pubchem_id = props[1].text
                break # romper una vez que encuentre el pubchem 
        dbids.append(drugbank_ID)
        pubchemids.append(pubchem_id)
    #
    dic_cid_dbid = dict((zip(pubchemids, dbids)))
    # first element is nan, delete
    dic_cid_dbid[list(dic_cid_dbid.keys())[0]]
    dic_cid_dbid.pop(list(dic_cid_dbid.keys())[0])
    return dic_cid_dbid

def get_drug_drug_coordinates(drug_entry, drug_drug_coodinates):
    drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
    all_interactions =  drug_entry.findall('.//{http://www.drugbank.ca}drug-interactions')[0]
    for inter in all_interactions: # aqui iteramos en targets
        coordinate = (drugbank_ID, list(inter)[0].text)	
        drug_drug_coodinates.append(coordinate)
    return drug_drug_coodinates

def select_drugbankid(x):
    '''
    This function returns the pubchem id depending on what code is available in data
    in order to retrieve the max possible information.
    '''
    if x['Pubchem_map_flat'] == x['Pubchem_map_stereo']:
        return x['Pubchem_map_flat']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and (x['Pubchem_map_stereo']) != np.nan ) :
        return x['Pubchem_map_stereo']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and (np.isnan(x['Pubchem_map_flat']) == False)):
        return x['Pubchem_map_flat']
    else:
        return np.nan

def read_fasta(path):
    names=[]
    seqs = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                names.append(line.strip().replace('>', ''))
            else:
                seqs.append(line.strip())
    return zip(names, seqs)

def write_fasta(path, target, seq):
    fl_name = os.path.join(path, target.replace(':', '_')+'.fasta')
    if os.path.exists(fl_name):
        logging.debug(f'File {fl_name} already exists')
        return fl_name
    with open(fl_name, 'w') as f:
        _ = f.write('>'+target+'\n'+seq+'\n')
    return fl_name

def create_remove_tmp_folder(path):
    if not os.path.exists(path):
        logging.info('Creating tmp folder: {}'.format(path))
        os.makedirs(path)
        return path
    else: 
        return path


def check_and_create_fasta(target, seq):
    global PATH
    fasta1 = os.path.join(PATH, target.replace(':', '_')+'.fasta')
    if not os.path.exists(fasta1):
        fasta1 = write_fasta(PATH, target, seq)
    return fasta1


def get_amino_uniprot(proteinID):
    r = requests.get(f'https://www.uniprot.org/uniprot/{proteinID}.fasta')
    if r.status_code == 200 and r.text:
        return (proteinID, ''.join(r.text.split('\n')[1:]))
    else: 
        #print('Protein sequence not found in uniprot database')
        return (proteinID, None)


def get_pairwise_tanimoto(smiles1,smiles2, dic): #dic == smile2fp
    try:
        for smile in [smiles1, smiles2]:
            if not smile in dic:
                mol1 = Chem.MolFromSmiles(str(smile))
                fp1  = Chem.RDKFingerprint(mol1)
                dic[smile] = fp1
        tani = DataStructs.FingerprintSimilarity(dic[smiles1],dic[smiles2]) #pairwise similarity
        return tani
    except:
        return None


#
def check_drug(drug_entry):
    try:
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        kegg_index = [ _.text for _ in drug_entry.findall('.//{http://www.drugbank.ca}resource')].index('KEGG Drug')
        kegg_ID = drug_entry.findall('.//{http://www.drugbank.ca}identifier')[kegg_index].text
        return drugbank_ID, kegg_ID
    except ValueError:
        return None, None

def get_dict_kegg2db():
    logging.debug('Reading DrugBank xml file...')
    tree = ET.parse('../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    logging.debug('Succesfully read!')
    root = tree.getroot()
    dbids = []
    keggids = []
    for drug_entry in tqdm(root):
        drugbank_ID, kegg_ID = check_drug(drug_entry)
        if drugbank_ID and kegg_ID:
            dbids.append(drugbank_ID)
            keggids.append(kegg_ID)
    kegg2db = dict(zip(keggids, dbids))
    return kegg2db

def get_smiles(drug_entry):
    '''
    Get list of drugs with smiles from DrugBank xml, 
    for SMILES that are not in DrugBank retrieves them from PubChem
    Then, check if we can create a fingerprint (if not, do not include)
    '''
    drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
    smiles = None
    fp = None
    for props in drug_entry.findall('.//{http://www.drugbank.ca}property'):
        for prop in props: 
            if(prop.text == 'SMILES'):
                smiles = props[1].text
                break
    if not smiles:
        for exids in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
            for ids in exids:
                if(ids.text == 'PubChem Compound'): 
                    pubchem_id = exids[1].text
                    smiles = get_compound_pubchem(pubchem_id)
                    break
    if not smiles:
        return(drugbank_ID, None)
    elif Chem.MolFromSmiles(str(smiles)):
        return(drugbank_ID, smiles)
    return(drugbank_ID, None)


def get_dict_hsa2uni():
    '''
    http://rest.kegg.jp/conv/hsa/uniprot 
    return a dict to change from hsa to uniprot ID
    '''
    r = requests.get('http://rest.kegg.jp/conv/hsa/uniprot')
    rlines = r.text.split('\n')
    up, hsa = [], []
    for line in rlines:
        if line != '':
            upi, hsai = line.split('\t')
            up.append(upi.lstrip('up:'))
            hsa.append(hsai)
    return dict(zip(hsa, up))


def get_dict_kegg2pubchemsid():
    '''
    http://rest.kegg.jp/conv/drug/pubchem 
    return a dict to change from kegg to PubChem ID SID!!
    '''
    r = requests.get('http://rest.kegg.jp/conv/drug/pubchem')
    rlines = r.text.split('\n')
    pubs, keggs = [], []
    for line in rlines:
        if line != '':
            pub, kegg = line.split('\t')
            pubs.append(pub.lstrip('pubchem:'))
            keggs.append(kegg.lstrip('dr:'))
    return dict(zip(keggs, pubs))

# new from get_SW_score_Yamanishi.py in DTI2Vec
def get_SW_score(pair1, pair2, tmp_path):
    target1, _ = pair1
    target2, _ = pair2
    fasta1 = os.path.join(tmp_path, target1.replace(':', '_')+'.fasta')
    fasta2 = os.path.join(tmp_path, target2.replace(':', '_')+'.fasta')
    args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
            '-asequence', fasta1 , '-bsequence', fasta2, 
            '-gapopen', '10.0', '-gapext', '0.5', 
            '-stdout']
    try:
        score = sp.Popen(args, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.DEVNULL)
        score = score.communicate(b"\n")
        score = score[0].decode().split('\n')
        score = extract_score(score)
        return score
    except:
        logging.warning(f'Not able to compute SW score for : {target1}, {target2}')

def extract_score(score):
    score = [line for line in score if line.startswith('# Score:')]
    if score:
        return float(score[0].split()[-1])
    else:
        return None

def write_all_fastas(fastas, path):
    for header, seq in fastas:
        write_fasta(path, header, seq)
    logging.info('All fastas written')



def drugname_drugbankid():
    '''
    This function creates a dictionary from drugbank name to DrugbankID
    (parsin the DrugBank xml file)
    '''
    logging.debug("Reading xml")
    tree = ET.parse('../../DB/Data/DrugBank/full_database.xml')
    root = tree.getroot()
    logging.debug("creating dict")
    drugbank_dic = {} # dict with generic or product name as key and drugbank_id as value
    for drug_entry in tqdm(root):
        drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
        #if drugbank_ID in all_entries:
        name = drug_entry.find('{http://www.drugbank.ca}name').text.lower()
        prod = drug_entry.find('{http://www.drugbank.ca}products')
        prod_names = set([brandname.find('{http://www.drugbank.ca}name').text.upper() for brandname in prod])
        if name:
            drugbank_dic[name] = drugbank_ID
        if len(prod_names) >= 1:
            for prod in prod_names:
                drugbank_dic[prod] = drugbank_ID
    return drugbank_dic



def get_dic_cid2sid(drug_list, chunk_size = 100, time_sleep = 1):
    '''
    *** explain ***
    '''
    n = chunk_size
    list_drug_list = [drug_list[i:i + n] for i in range(0, len(drug_list), n)]
    d = {}
    for lts in tqdm(list_drug_list):
        pclist = ','.join(lts)
        r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pclist}/sids/JSON')
        data = r.json()['InformationList']['Information'] #data is a list of dictionaries
        for dicti in data:
            if dicti.get('SID'):
                cid = str(dicti.get('CID'))
                sids = [str(si) for si in dicti.get('SID')]
                d[cid] = sids
        time.sleep(time_sleep)
    return d


def get_dic_cid2syn(drug_list, chunk_size = 100, time_sleep = 1):
    '''
    *** explain ***
    '''
    n = chunk_size
    list_drug_list = [drug_list[i:i + n] for i in range(0, len(drug_list), n)]
    d = {}
    for lts in tqdm(list_drug_list):
        pclist = ','.join(lts)
        r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pclist}/synonyms/JSON')
        data = r.json()['InformationList']['Information'] #data is a list of dictionaries
        for dicti in data:
            if dicti.get('Synonym'):
                cid = str(dicti.get('CID'))
                syns = [str(si) for si in dicti.get('Synonym')]
                d[cid] = syns
        time.sleep(time_sleep)
    return d




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





def get_dict_dtis_yamanishi():
    GEN_PATH = '../../DB/Data/'
    PATH_E = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt')
    PATH_NR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt')
    PATH_GPCR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt')
    PATH_IC = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt')
    list_yam = [PATH_E, PATH_NR, PATH_GPCR, PATH_IC]
    names_yam = ['E', 'NR', 'GPCR', 'IC']
    dics_yam = {}
    kegg2sid = get_Kegg_pubchem_dict()
    hsa2uni = get_dict_hsa2uni()
    for yam_db,name in zip(list_yam,names_yam):
        colnames_data = ['Kegg_ID', 'Gene']
        df = pd.read_csv(yam_db, header = None, names = colnames_data, index_col=False, sep='\t')
        df['sid'] = df.Kegg_ID.map(kegg2sid)
        sid2cid = get_cid_from_sid(df['sid'].dropna().unique().tolist())
        df['PubChemID'] = df.sid.map(sid2cid)
        # proteins 
        df['Gene'] = df['Gene'].map(lambda x: x[:3] + ":" + x[3:])
        df['UniprotID'] = df['Gene'].map(hsa2uni)
        df = df.drop(columns = ['Kegg_ID', 'Gene', 'sid'])
        df = df.dropna().drop_duplicates()
        df['PubChemID'] = df.PubChemID.astype(int).astype(str)
        dics_yam[name] = df
    return dics_yam


def get_dtis_drugbank():
    GEN_PATH = '../../DB/Data/'
    PATH_DRUGBANK = os.path.join(GEN_PATH, 'DrugBank/DrugBank_DTIs.tsv')
    df = pd.read_csv(PATH_DRUGBANK, sep='\t') 
    df.columns = ['DrugBankID', 'UniprotID']
    df = df.drop_duplicates()
    dic_cid_dbid = pubchem_to_drugbankid()
    dic_dbid_cid = {v: k for k, v in dic_cid_dbid.items()}
    df['PubChemID'] = df.DrugBankID.map(dic_dbid_cid)
    df = df.dropna()
    df = df.drop_duplicates()
    df = df.drop(columns = 'DrugBankID')
    df = df[['PubChemID', 'UniprotID']]
    return df


def get_dtis_davis():
    GEN_PATH = '../../DB/Data/'
    PATH_DAVIS = os.path.join(GEN_PATH, 'Davis_et_al/tdc_package_preprocessing/DAVIS_et_al_w_labels.tsv')
    davis = pd.read_csv(PATH_DAVIS, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Label'])
    davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
    mod_davis = davis.copy()
    # then drop those with 0
    mod_davis = mod_davis.drop(mod_davis[mod_davis.Label == 0].index)
    mod_davis = mod_davis.drop(columns=['SMILES','Label', 'Sequence'])
    ## PROTEIN IDENTIFIERS
    data_path_biomart = os.path.join(GEN_PATH, 'cross_side_information_DB/bioMART/mart_export_expanded.txt')
    biomart = pd.read_csv(data_path_biomart, sep='\t', usecols=['UniProtKB Gene Name symbol', 'UniProtKB Gene Name ID'])
    biomart = biomart.dropna().drop_duplicates()
    genename2geneid = dict(zip(biomart['UniProtKB Gene Name symbol'].tolist(), biomart['UniProtKB Gene Name ID'].tolist()))
    mod_davis['UniprotID'] = mod_davis['GeneName'].map(genename2geneid)
    mod_davis = mod_davis.dropna().drop_duplicates()
    df = mod_davis[['PubChemID', 'UniprotID']]
    return df


def get_bindingdb_dtis():
    db_file_path = '../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
    bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Y'])
    bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID', 'Target Sequence': 'Sequence'}, axis=1)
    mod_bind = bindingdb.copy()
    threshold = 30
    mod_bind['Label'] = [1 if x < threshold else 0 for x in mod_bind['Y']]
    mod_bind = mod_bind.drop(mod_bind[mod_bind.Label == 0].index)
    mod_bind = mod_bind.drop(columns=['Y', 'Label', 'Sequence', 'SMILES'])
    mod_bind['PubChemID'] = mod_bind.PubChemID.astype(int).astype(str)
    df = mod_bind
    df = df.dropna().drop_duplicates()
    return df


def get_biosnap_dtis():
    dti_file_path = '../../DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
    dti = pd.read_csv(os.path.join(os.getcwd(), dti_file_path), sep='\t', comment='#', header=None) 
    dti.columns = ['DrugBankID', 'UniprotID']
    # get dict and inverse it
    dic_cid_dbid = pubchem_to_drugbankid()
    dic_dbid_cid = {v: k for k, v in dic_cid_dbid.items()}
    dti['PubChemID'] = dti.DrugBankID.map(dic_dbid_cid)
    df = dti[['PubChemID', 'UniprotID']].dropna().drop_duplicates()
    return df



### get dtis original xx
class original_dtis:
        def __init__(self):
                self.GEN_PATH = '../../DB/Data/'
        #
        def biosnap(self):
                GEN_PATH = self.GEN_PATH
                PATH_BIOSNAP = 'BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
                dti_file_path = os.path.join(GEN_PATH, PATH_BIOSNAP)
                dti = pd.read_csv(os.path.join(os.getcwd(), dti_file_path), sep='\t', comment='#', header=None) 
                dti.columns = ['DrugBankID', 'UniprotID']
                dti = dti.dropna().drop_duplicates()
                return dti
        #
        def bindingdb(self):
                GEN_PATH = self.GEN_PATH
                PATH_BINDING = 'BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
                db_file_path = os.path.join(GEN_PATH, PATH_BINDING)
                bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Y'])
                bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID', 'Target Sequence': 'Sequence'}, axis=1)
                mod_bind = bindingdb.copy()
                threshold = 30
                mod_bind['Label'] = [1 if x < threshold else 0 for x in mod_bind['Y']]
                mod_bind = mod_bind.drop(mod_bind[mod_bind.Label == 0].index)
                mod_bind = mod_bind.drop(columns=['Y', 'Label', 'Sequence', 'SMILES'])
                mod_bind['PubChemID'] = mod_bind.PubChemID.astype(int).astype(str)
                df = mod_bind
                df = df.dropna().drop_duplicates()
                return df
        #
        def davis(self):
                GEN_PATH = self.GEN_PATH
                PATH_DAVIS = os.path.join(GEN_PATH, 'Davis_et_al/tdc_package_preprocessing/DAVIS_et_al_w_labels.tsv')
                davis = pd.read_csv(PATH_DAVIS, sep="\t", header=0, usecols=['Drug_ID', 'SMILES', 'Target_ID', 'Target Sequence', 'Label'])
                davis = davis.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'GeneName', 'Target Sequence': 'Sequence'}, axis=1)
                mod_davis = davis.copy()
                # then drop those with 0
                mod_davis = mod_davis.drop(mod_davis[mod_davis.Label == 0].index)
                mod_davis = mod_davis.drop(columns=['SMILES','Label', 'Sequence'])
                mod_davis = mod_davis.dropna().drop_duplicates()
                return mod_davis
        #
        def drugbank(self):
                GEN_PATH = self.GEN_PATH
                PATH_DRUGBANK = os.path.join(GEN_PATH, 'DrugBank/DrugBank_DTIs.tsv')
                df = pd.read_csv(PATH_DRUGBANK, sep='\t') 
                df.columns = ['DrugBankID', 'UniprotID']
                df = df.dropna().drop_duplicates()
                return df
        #
        def dict_yamanishi(self):
                GEN_PATH = self.GEN_PATH
                PATH_E = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt')
                PATH_NR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt')
                PATH_GPCR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt')
                PATH_IC = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt')
                list_yam = [PATH_E, PATH_NR, PATH_GPCR, PATH_IC]
                names_yam = ['E', 'NR', 'GPCR', 'IC']
                dics_yam = {}
                for yam_db,name in zip(list_yam,names_yam):
                        colnames_data = ['Kegg_ID', 'Gene']
                        df = pd.read_csv(yam_db, header = None, names = colnames_data, index_col=False, sep='\t')
                        df = df.dropna().drop_duplicates()
                        dics_yam[name] = df
                return dics_yam

