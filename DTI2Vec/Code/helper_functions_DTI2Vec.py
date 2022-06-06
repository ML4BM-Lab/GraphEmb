import requests
import os, sys, uuid, re
import pandas as pd
import subprocess as sp
import multiprocessing as mp
from itertools import repeat
from re import search
from scipy.fft import idct
from tqdm import tqdm
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler
import logging
import xml.etree.ElementTree as ET
from pubchempy import Compound

def get_DB_name(path):
    """
    This function returns the name of the DB.
    """
    DB_NAMES = ['BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'DrugBank' ,'E', 'GPCR', 'IC', 'NR']
    for db in DB_NAMES:
        if search(db, path):
            logging.info(f'Database: {db}')
            if db in ['E', 'GPCR', 'IC', 'NR']:
                db = os.path.join('Yamanashi_et_al_GoldStandard', db)
                if not os.path.exists(os.path.join('./../Data', db)):
                    logging.info(f'Creating direcotry: ./../Data/{db}')
                    os.mkdir(os.path.join('./../Data', db))
                return db
            else:
                if not os.path.exists(os.path.join('./../Data', db)):
                    logging.info(f'Creating direcotry: ./../Data/{db}')
                    os.mkdir(os.path.join('./../Data', db))
                return db
    logging.error(f'Database: {db} not found')
    sys.exit('Please provide a valid database')

def read_dtis_Yamanashi(path):
    """"
    Read the dti file and return the dti along with 
    the numerical dictionary sorted by name
    """
    with open(path, 'r') as f:
        dtis = f.readlines()
    dtis = [entry.strip().split('\t') for entry in dtis]
    all_elements = list(set(element for entry in dtis for element in entry))
    logging.info(f'Number of unique nodes: {len(all_elements)}')
    logging.info(f'Number of unique edges: {len(dtis)}')
    dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
    return dtis, dictionary 

def read_and_extract_DAVIS(path):
    drugs = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('\t'):
                drug = line.split('\t')[1]
                drugs.append(drug)
    return drugs

def read_dtis_DrugBank(path):
    """"
    Read the dti file and return the dti along with 
    the numerical dictionary sorted by name
    """
    with open(path, 'r') as f:
        _ = next(f)
        dtis = f.readlines()
    dtis = [entry.strip().split('\t') for entry in dtis]
    all_elements = list(set(element for entry in dtis for element in entry))
    logging.info(f'Number of unique nodes: {len(all_elements)}')
    logging.info(f'Number of unique edges: {len(dtis)}')
    dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
    return dtis, dictionary 

def read_dtis_BIOSNAP(path):
    """"
    Read the dti file and return the dti along with 
    the numerical dictionary sorted by name
    """
    with open(path, 'r') as f:
        _ = next(f)
        dtis = f.readlines()
    dtis = [entry.strip().split('\t') for entry in dtis]
    all_elements = list(set(element for entry in dtis for element in entry))
    logging.info(f'Number of unique nodes: {len(all_elements)}')
    logging.info(f'Number of unique edges: {len(dtis)}')
    dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
    return dtis, dictionary 

def get_seqs_DAVIS(path):
    """
    This function reads the database and returns the targets 
    """
    targets = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('\t'):
                line = line.split('\t')
                targets.append((line[3], line[4]))
        return targets

def read_dtis_DAVIS(path):
    """"
    Read the dti file and return the dti along with 
    the numerical dictionary sorted by name
    """
    with open(path, 'r') as f:
        _ = next(f)
        dtis = f.readlines()
    dtis = [(entry.split('\t')[1], entry.split('\t')[3])  for entry in dtis]
    all_elements = list(set(element for entry in dtis for element in entry))
    logging.info(f'Number of unique nodes: {len(all_elements)}')
    logging.info(f'Number of unique edges: {len(dtis)}')
    dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
    return dtis, dictionary

def read_dtis_Binding(path):
    """"
    Read the dti file and return the dti along with 
    the numerical dictionary sorted by name
    """
    with open(path, 'r') as f:
        _ = next(f)
        dtis = f.readlines()
    dtis = [(entry.split('\t')[1], entry.split('\t')[3])  for entry in dtis]
    all_elements = list(set(element for entry in dtis for element in entry))
    logging.info(f'Number of unique nodes: {len(all_elements)}')
    logging.info(f'Number of unique edges: {len(dtis)}')
    dictionary = {node : str(index)  for index, node  in enumerate(sorted(all_elements))}
    return dtis, dictionary

def get_primary_drugbank_ID(drugbank_ID, drug_entry):
    drugbank_IDs = drug_entry.findall('{http://www.drugbank.ca}drugbank-id')
    attribs = [elem.attrib for elem in drugbank_IDs]
    primary = [attrib for attrib in attribs if attrib][0]
    if primary['primary'] == 'true':
        drugbank_ID = drugbank_IDs[attribs.index(primary)].text
    return drugbank_ID

def get_drugs_smiles_from_DrugBank(drugs):
    logging.info('Loading DrugBank...')
    tree = ET.parse('./../../DB/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
    root = tree.getroot()
    drugs_smiles = []
    for drug_entry in tqdm(root):
        drugbank_IDs = drug_entry.findall('{http://www.drugbank.ca}drugbank-id')
        drugbank_IDs = [name.text for name in drugbank_IDs]
        drugbank_ID = [real_name for real_name in drugbank_IDs if real_name in drugs]
        if drugbank_ID:
            if len(drugbank_ID) > 1:
                drugbank_ID = get_primary_drugbank_ID(drugbank_ID, drug_entry)
            else:
                drugbank_ID = drugbank_ID[0]
                props = drug_entry.findall('.//{http://www.drugbank.ca}property')
                smile = [prop for prop in props if prop[0].text == 'SMILES']
                if smile:
                    # if smile in drugbank data
                    smile = smile[0][1].text
                else:
                    # if smile not in drugbank data retrieve it from chembl
                    external_ids = drug_entry.findall('.//{http://www.drugbank.ca}external-identifier')
                    pubchem_id = [external_id for external_id in external_ids if external_id[0].text == 'PubChem Compound']
                    if pubchem_id:
                        pubchem_id = pubchem_id[0][1].text
                        smile = Compound.from_cid(pubchem_id).canonical_smiles
                drugs_smiles.append((drugbank_ID, smile))
    return drugs_smiles

def get_seqs_BindingDB(path):
    """
    This function reads the database and returns the targets 
    """
    targets = []
    with open(path, 'r') as f:
        _ = next(f)
        for line in f:
            line = line.split('\t')
            targets.append((line[3], line[4]))
        return targets

def read_and_extract_BINDING_smiles(path):
    """
    Read the BINDING file and return the smiles
    """
    with open(path, 'r') as f:
        _ = next(f)
        smiles = f.readlines()
    smiles = [entry.strip().split('\t') for entry in smiles]
    smiles = [(entry[1], entry[2]) for entry in smiles]
    smiles = list(set(smiles))
    return smiles

def get_BIOSNAP_drugs(path):
    with open(path, 'r') as f:
        _ = next(f)
        drugs = f.readlines()
    drugs = [entry.strip().split('\t') for entry in drugs]
    drugs = [entry[0] for entry in drugs]
    drugs = list(set(drugs))
    return drugs	



def write_dtis(dtis, path):
    result_ID = str(uuid.uuid4())
    result_file = os.path.join(path, result_ID+'_coded_dti.txt')
    TAB = '\t'
    NL = '\n'
    dtis = [TAB.join(entry) for entry in dtis]
    with open(result_file, 'w') as outfl:
        outfl.writelines(NL.join(dtis))
    logging.debug(f'Coded DTIs written at {result_file}') 
    return result_file

def write_dict(dictionary, file_path):
    with open(file_path, 'w') as f:
        for key, value in dictionary.items():
            f.write(f'{key}\t{value}\n')
    logging.debug(f'Coded DTIs written at {file_path}') 

def write_embedddings_with_name(embeddings, node_dict):
    with open(embeddings , 'r') as infl :
        stats = next(infl).strip().split()
        embs = [line.strip().split() for line in infl]

    assert len(embs)    == int(stats[0]), 'Number of embeddings does not match'
    assert len(embs[0]) == (int(stats[1]) +1 ), 'Number of dimensions does not match'
    inv_map = {v: k for k, v in node_dict.items()}
    SEP = ' '
    named_nodes = embeddings.replace('.txt', '_with_name.txt')
    with open(named_nodes, 'w') as outfl:
        for emb in embs:
            node = emb[0]
            emb = SEP.join(emb[1:])
            outfl.write(f'{inv_map[node]}{SEP}{emb}\n')
    logging.debug(f'Embeddings written at {named_nodes}')

def write_smile_per_file(drug_smiles, db_name):
    SCRATCH_DISK = '/media/scratch_ssd/tmp/'
    if os.path.isdir(SCRATCH_DISK):
        tmp_path = create_remove_tmp_folder(os.path.join(SCRATCH_DISK , db_name))
    else:
        tmp_path = create_remove_tmp_folder(os.path.join('/tmp/SIMCOMP_alternative' , db_name))
    for id, smile in tqdm(drug_smiles, desc='Writing smiles'):
        with open(os.path.join(tmp_path, f'{id}.smi'), 'w') as f:
            _ = f.write(smile)
    return tmp_path

def create_remove_tmp_folder(path):
    if not os.path.exists(path):
        logging.info('Creating tmp folder: {}'.format(path))
        os.makedirs(path)
        return path
    else: 
        return path

def read_and_extract_targets(path):
    """
    This function reads the database and returns the targets.
    """
    targets = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                _, target = line.split('\t')
                targets.append(target.strip())
    return targets

def getamino_KEGG(protein):
    r = requests.get(f'http://rest.kegg.jp/get/{protein}/aaseq')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

def replace_and_get_AA(ID):
    """
    This function replaces the 'hsa' in the identifier with 'hsa:'.
    """
    ID = ID.replace('hsa', 'hsa:')
    return (ID, getamino_KEGG(ID))

def get_drug_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

def get_yamanashi_subDB(path):
    return re.search(r'(?<=\/)[\w]+', path).group()

def retrieve_sequences_from_UniProt(ID):
    """
    This function retrieves the AA sequence from uniprot.
    """
    r = requests.get(f'https://www.uniprot.org/uniprot/{ID}.fasta')
    if r.status_code == 200:
        aminoseq = ''.join(r.text.split('\n')[1:])
    else:
        logging.error(f'{ID} not found')
        return (ID, None)
    return (ID, aminoseq)

def extract_score(results_file):
    with open(results_file, 'r') as f:
        for line in f:
            if not line.startswith('# Score:'):
                continue
            else:
                return float(line.split()[-1])

def write_fasta(path, target, seq):
    fl_name = os.path.join(path, target.replace(':', '_')+'.fasta')
    if os.path.exists(fl_name):
        logging.debug(f'File {fl_name} already exists')
        return fl_name
    with open(fl_name, 'w') as f:
        _ = f.write('>'+target+'\n'+seq+'\n')
    return fl_name

def check_and_create_fasta(target, seq):
    global PATH
    fasta1 = os.path.join(PATH, target.replace(':', '_')+'.fasta')
    if not os.path.exists(fasta1):
        fasta1 = write_fasta(PATH, target, seq)
    return fasta1

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

def check_and_create_folder(db_name):
    if not os.path.exists(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/', db_name)):
        os.mkdir(os.path.join('/home/margaret/data/jfuente/DTI/Input4Models/DTI2Vec/', db_name))

def read_AA_sequences(path):
    # read a fasta file and return a dictionary with the target as key and the sequence as value
    AA_sequences = {}
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                target = line.strip().replace('>', '')
            else:
                AA_sequences[target] = line.strip()
    return list(AA_sequences.items())

def write_all_fastas(fastas, path):
    for header, seq in fastas:
        write_fasta(path, header, seq)
    logging.info('All fastas written')

def get_dti(path):
    with open(path, 'r') as f:
        dtis = f.readlines()
        dtis = [line.strip().split('\t') for line in dtis]
    return dtis

def get_admat_from_dti(edges):
    # get the nodes:
    horizontal_nodes = list(set([edge[0] for edge in edges]))
    horizontal_nodes.sort()
    horizontal_node_position = {node: i for i, node in enumerate(horizontal_nodes)}
    vertical_nodes = list(set([edge[1] for edge in edges]))
    vertical_nodes.sort()
    vertical_node_position = {node: i for i, node in enumerate(vertical_nodes)}
    # get the matrix dimensions
    horizontal_length = len(horizontal_nodes)
    vertical_length = len(vertical_nodes)
    adjacency = [[0]*horizontal_length for _ in range(vertical_length)]
    for orig, dest in edges:
        adjacency[vertical_node_position.get(dest)][horizontal_node_position.get(orig)] = 1
    adjacency  = pd.DataFrame(adjacency, columns=horizontal_nodes, index=vertical_nodes)
    return adjacency

def write_edges(edges, path):
    # with open(path, 'w') as f:
    # 	for edge in edges:
    # 		f.write('\t'.join(edge)+'\n')
    # we need the edges to be in the format: target -- drug
    with open(path, 'w') as f:
        for edge1, edge2 in edges:
            f.write(edge2 +'\t'+ edge1+'\n')

def get_SMILES_from_Pubchem_web_batch(drugs):
    smiles = []
    # split the list in chunks of 100 to make requests
    size = 100
    drugs = [str(drug) for drug in drugs]
    split_list = lambda big_list, x: [
        big_list[i : i + x] for i in range(0, len(big_list), x)
    ]
    drug_chunks = split_list(drugs, size)
    for chunk in tqdm(drug_chunks):
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
        if smile:
            smiles.append((cid, smile[0]))
    return smiles

def write_smiles(db_name, drugs_smiles):
    with open(os.path.join('/./../Data/DTI2Vec/', db_name, 'drugs_smiles.tsv'), 'w') as f:
        for drug_id, smiles in drugs_smiles:
            if smiles:
                _ = f.write(f'{drug_id}\t{smiles}\n')
