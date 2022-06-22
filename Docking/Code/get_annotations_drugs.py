import os
import pandas as pd
import logging
from tqdm import tqdm
import helper_functions as hf
import requests
from tqdm import tqdm
import numpy as np
import requests
import subprocess as sp
from itertools import repeat
from rdkit import Chem

#### Return all drugs in PubChemID

# Calculate similarity with tanimoto
# fingerprint by default
# Get Family anotations from ClassyFire


###########

def get_SMILES_from_Pubchem_web_batch(drugs):
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
            if smile:
                try:
                    mol1 = Chem.MolFromSmiles(str(smile[0]))
                    fp1  = Chem.RDKFingerprint(mol1)
                    smiles.append((cid, smile[0]))
                except:
                    logging.info(f'Error for pubchemid {cid}')
    return smiles




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


###########

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)

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




## yamanishi

GEN_PATH = '../../DB/Data/'
PATH_E = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt')
PATH_NR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt')
PATH_GPCR = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt')
PATH_IC = os.path.join(GEN_PATH, 'Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt')
list_yam = [PATH_E, PATH_NR, PATH_GPCR, PATH_IC]
#for yam_db in list_yam:

yam_db = PATH_E
print(yam_db)
# yam
colnames_data = ['Kegg_ID', 'Gene']
df = pd.read_csv(yam_db, header = None, names = colnames_data, index_col=False, sep='\t')
df

kegg2db = hf.get_dict_kegg2db()
df['DrugBank_ID'] = df['Kegg_ID'].map(kegg2db)

dic_cid_dbid =  hf.pubchem_to_drugbankid()
dic_dbid_cid = {v: k for k, v in dic_cid_dbid.items()}

drugs = df['DrugBank_ID'].unique().tolist()

drugs = df['Kegg_ID'].unique().tolist()

drugs_pc = [dic_dbid_cid.get(i, None) for i in drugs] # 
drugs_pc = list(filter(None, drugs_pc)) 



import requests


url = 'https://cts.fiehnlab.ucdavis.edu/rest/convert/KEGG/PubChem%20CID/D03752'

for drug in drugs:
    drug = str(drug)
    resp = requests.get(f'https://cts.fiehnlab.ucdavis.edu/rest/convert/KEGG/PubChem%20CID/{drug}')
    if resp.status_code == 200:
        resp.json()

# https://rest.kegg.jp/conv/drug/pubchem
# curl https://cts.fiehnlab.ucdavis.edu/rest/convert/KEGG/PubChem%20CID/D03752

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

import re
dictff = get_Kegg_pubchem_dict()


np.savetxt('kegg_ids.txt', list(dictff.keys()), fmt='%s')

###################
# Load Drugs & SMILES
drugs_drugbank = get_drugs_drugbank() # only drugbank
drugs_davis = get_drugs_Davis()
drugs_binding = get_drugs_binding()
drugs_biosnap = get_drugs_biosnap()

all_drugs = list(set(drugs_drugbank + drugs_davis + drugs_binding + drugs_biosnap))

logging.info(f'Total number of drugs: {len(all_drugs)}')

list_drug_smiles_inch = get_SMILES_n_InCh_from_Pubchem_web_batch(all_drugs) 


##### SMILES

############### TANIMOTO SCORE FOR CHEMICAL SIMILITUDE

list_drugs = [code for code, _, _ in list_drug_smiles_inch]
list_smiles =  [smiles for _, smiles,_ in list_drug_smiles_inch] #[dict_drugid_smiles[i] for i in list_drugs] # 
list_inchi =  [inchi for _,_,inchi in list_drug_smiles_inch]
#dict_drugid_smiles = dict(zip(list_drugs, list_smiles))

any(elem is None for elem in list_inchi)


###
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

#
df_all_Sim_Tani = pd.DataFrame(all_Sim_Tani, columns= list_drugs, index = list_drugs) 

logging.debug(f'Drug Similarity Matrix Shape {df_all_Sim_Tani.shape}')
# save pickle
# save csv for model

###########################
## Annotation


############################
###### REPRS
import matplotlib.pyplot as plt
import seaborn as sns

df = df_all_Sim_Tani.astype(float)
# df.to_pickle('test_tani.pkl')
# df = pd.read_pickle('test_tani.pkl')

plt.clf()
plt.title('Tanimoto Similitude (DrugBank)')
sns.heatmap(df, xticklabels=False, yticklabels=False)
#sns.clustermap(df, metric='euclidian')
plt.savefig('../Results/drug_test_2.png',dpi=300)

# plt.imshow(np.array(df.values.tolist()).astype('float'))
# plt.imshow(df, cmap=plt.cm.get_cmap("Reds"), interpolation="nearest")
# plt.colorbar()

plt.clf()
plt.title('Tanimoto Similitude Cluster (DrugBank)')
#sns.heatmap(df, xticklabels=False, yticklabels=False)
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)

sns.clustermap(df, metric='euclidean',  cmap = pal, xticklabels=False, yticklabels=False)

plt.savefig('../Results/drug_test_2_clust_nolab.png',dpi=400)



#####


# df.index.tolist()[0]
# means = []
# for col in df.index.tolist():
#     means.append(df[col].drop(col).mean())

# plt.clf()
# plt.title('hist drugbank mean')
# #plt.hist(means, bins=100)
# #sns.displot(means, kind='kde', cut=0)
# sns.displot(means, kde=True)

plt.clf()
plt.title('hist drugbank')
keep = np.triu(np.ones(df.shape)).astype('bool').reshape(df.size)
vals = df.stack()[keep].to_numpy()
#sns.displot(vals, kind='kde')
sns.histplot(vals, kde=True)
plt.savefig('../Results/drug_test_hist.png',dpi=300)

# heatmap of full dataframe



#### Classify for annotation
# http://classyfire.wishartlab.com/queries/new



#### get annotations
annotations = []
retrieve = ['kingdom', 'superclass', 'class', 'subclass']
for i in tqdm(range(len(list_inchi)), desc='Annotating drugs'):
    elmt = list_inchi[i]
    url = f'http://classyfire.wishartlab.com/entities/{elmt}.json'
    response = requests.get(url)
    results_ = []
    if response.status_code == 200:
        json = response.json()
        #json.keys()
        for item in retrieve:
            t = json.get(item, None)
            if t: 
                result = t.get('name', None)
                results_.append(result)
    results_.insert(0, list_drugs[i])
    annotations.append(tuple(results_))



df_annot = pd.DataFrame.from_records(annotations, columns=['PubChemID', 'kingdom', 'superclass', 'class', 'subclass'])
df_annot=df_annot.fillna('-')
#########################
##### test clustermap
from matplotlib.patches import Patch

plt.clf()
#sns.heatmap(df, xticklabels=False, yticklabels=False)
#pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)
#sns.clustermap(df, metric='euclidean',  cmap = pal, xticklabels=False, yticklabels=False)
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)

#dict_classes = {superclass_unique[i]:i+1 for i in range(len(df_annot.superclass.unique()))}

types = matrix['superclasse']
colors_type = sns.color_palette("ch:s=.70,rot=-.70", 8)
sorted_types = types.unique().tolist()
sorted_types.sort()
lut_type = dict(zip(sorted_types,colors_type ))
row_colors = types.map(lut_type)
# labels = np.random.random_integers(0,5, size=50)
# lut = dict(zip(set(labels), sns.hls_palette(len(set(labels)), l=0.5, s=0.8)))
# row_colors = pd.DataFrame(labels)[0].map(lut)

#sns.clustermap(matrix.iloc[:,:-1],yticklabels=False,xticklabels=False)
sns.clustermap(matrix.iloc[:,:-1], cmap = pal,row_colors=row_colors, yticklabels=False, xticklabels=False)
# handles = [Patch(facecolor=lut[name]) for name in lut]
# plt.legend(handles, lut, title='Types',
#            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
plt.savefig('../Results/drug_test_annot_2.png',dpi=400)


matrix = df
drug2superclass  = dict(zip(df_annot.PubChemID, df_annot.superclass))
matrix['superclasse'] = matrix.index.map(drug2superclass)
matrix