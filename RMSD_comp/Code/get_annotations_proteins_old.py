import requests
from requests.adapters import HTTPAdapter, Retry
import re
import logging
import numpy as np
import pandas as pd
import time
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
###

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)
###


def get_next_link(headers):
    # Code from Uniprot
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    # Code from Uniprot
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)



FILE_ALL_PROT = '../Data/pkls/annot_uniprot_all_raw.pkl'


if not os.path.isfile(FILE_ALL_PROT):

    ### for uniprot sesion
    logging.info('Retrieving annotations from Uniprot Rest API')
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    # load class
    class_list = np.loadtxt('../Data/keys.txt', dtype='str', delimiter='\n').tolist()
    class_list = [key.lower() for key in class_list]
    logging.debug("""keywords from: https://www.uniprot.org/keywords/KW-9992; 
                    deleting enzyme keywords because using EC Code for those""")

    # will retrieve a full pkl with all information for later only keeping the ones we want
    URL_HUMAN = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Clength%2Cec%2Ckeyword&format=tsv&query=%2A%20AND%20%28model_organism%3A9606%29&size=500'
    URL_WSTR = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Clength%2Cec%2Ckeyword&format=tsv&query=%2A%20AND%20%28proteins_with%3A1%29&size=500'

list_unis = []
for batch, total in get_batch(URL_WSTR):
    atime = time.time()
    for line in batch.text.splitlines()[1:]:
        uniprotid, length, ecode, kws = line.split('\t')
        ecode = ecode.split(';')[0] # only first
        kws = [x.lower() for x in kws.split(';')]
        kws = list(filter(lambda x: x != "", kws)) # check non empty strings
        if ecode:
            mol_funct = 'enzyme'
        elif any(key in class_list for key in kws):
            mol_funct = [word for word in class_list if word in kws][0]
        else:
            mol_funct = 'notAnnotated' 
        result = (uniprotid, int(length), ecode, mol_funct)
        list_unis.append(result)
    print(f'{len(list_unis)} / {total} : {(time.time()-atime):.4f} s')


    cols_all = ['UniprotID', 'length', 'EC', 'mol_func']
    df_all = pd.DataFrame.from_records(list_unis, columns=cols_all)
    df_all.to_pickle(FILE_ALL_PROT)

else:
    logging.info('Reading pkl')
    df_all = pd.read_pickle(FILE_ALL_PROT)



#list_proteins = np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()

# list of proteins w structure
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins_struct = list_pdbs + list_afold

#
# test difference
set(list_proteins_struct).difference(set(df_annot.UniprotID.unique().tolist()))
#

df_annot = df_all[df_all.UniprotID.isin(list_proteins_struct)].reset_index(inplace=False, drop=True)
df_annot.shape[0] == 4332

# removing those with freq < 10
frequencies = df_annot['mol_func'].value_counts()
condition = frequencies < 10
mask_obs = frequencies[condition].index
mask_dict = dict.fromkeys(mask_obs, 'other')
df_annot['mol_func'] = df_annot['mol_func'].replace(mask_dict) 

# Save distribution
plt.clf()
sns.set_theme(style="whitegrid")
ax = sns.countplot(x='mol_func', data=df_annot, order=df_annot.mol_func.value_counts().index)
ax.set_xticklabels(ax.get_xticklabels(),rotation = 35, size=7)
plt.savefig('../Results/test_distribution_molfunc.pdf')

## 
dic_alfa = dict(zip(list_afold, ['AFold'] * len(list_afold)))
dict_pdb = dict(zip(list_pdbs, ['PDB'] * len(list_pdbs)))
dic_source = {}
dic_source.update(dic_alfa)
dic_source.update(dict_pdb)

df_annot['source'] = df_annot.UniprotID.map(dic_source)


df_annot['EC'] = df_annot.EC[df_annot.EC != 0].apply(lambda x: x[:3])
df_annot['EC'] = df_annot.EC.replace('','-')

OUT_FILE = '../Data/pkls/df_annot_proteins_structure.pkl'
df_annot.to_pickle(OUT_FILE)



#### TEST  with curl
curl --form 'from=UniProtKB_AC-ID' \
     --form 'to=UniProtKB' \
     --form 'ids=B1AL88,X6RLY7' \
     https://rest.uniprot.org/idmapping/run



list_proteins = np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()

query = ",".join(list_proteins[:10])


idfrom = 'UniProtKB_AC-ID'
idto = 'UniProtKB'
ASKCOLUMNS = None

params = {
            'from': idfrom, 
            'to': 'ACC',
            'format': 'tab',
            'query': query,
            'columns':','.join(ASKCOLUMNS),
            }

r = requests.post('https://rest.uniprot.org/idmapping/run', data=params)

https://rest.uniprot.org/uniprotkb/accessions?accessions=A1L3X0%2CA6NGG8



