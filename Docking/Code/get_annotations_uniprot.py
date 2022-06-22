import os
import numpy as np
import pandas as pd
import glob
import logging
import requests
from tqdm import tqdm
import pandas as pd
import os
import logging
import re

#####
# 

logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)
##

##########################
########## CODE ##########
OUT_PATH = os.path.join(os.getcwd(), '../Data/pkls')

# Load objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]
list_proteins_struct = list_pdbs + list_afold

### instead # for all proteins
list_proteins = np.loadtxt('../Data/all_prot.txt', dtype='str').tolist()

class_list = np.loadtxt('keys.txt', dtype='str', delimiter='\n').tolist()
class_list = [key.lower() for key in class_list]
logging.debug("""keywords from: https://www.uniprot.org/keywords/KW-9992; 
                deleting enzyme keywords because using EC Code for those""")

logging.info(f'Retrieving information for {len(list_proteins)}')

# class_list = ['ion channel', 'g-protein coupled receptor', 'hormone',
#                 'cytokine', 'developmental protein', 'chaperone', 'activator'
#                 'antiviral protein', 'receptor', 'repressor','ribonucleoprotein', 
#                 'transport', 'transcription', 'transmembrane', 'membrane']

data_uniprot = pd.DataFrame(columns=['UniprotID', 'seq_len', 'Class', 'EC'])

for uniprotid in tqdm(list_proteins, desc='Looking up in Uniprot'):
    r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
    a = r.text.split('\n')
    SQ_info = None # sequence info
    check_SQ = any('SQ' in string for string in a)
    if check_SQ:
        SQ_info = [line.split() for line in a if 'SQ ' in line if check_SQ][0][2]
    ECode = None # enzyme code info
    check_E = any('EC=' in string for string in a)
    check_BR = any('BRENDA;' in string for string in a)
    if check_E: 
        E_Code_line = [line.split() for line in a if 'EC=' in line][0]
        ECode = [element for element in E_Code_line if '=' in element]
        ECode = ECode[0].replace('EC=', '').strip(';')
    elif check_BR:
        BR_Code_line = [line.split() for line in a if 'BRENDA;' in line]
        ECode = BR_Code_line[0][2].strip(';')
    check_KW = any('KW' in string for string in a) # find class
    if check_KW:
        KW_info = [line.replace('KW', '').replace('.', '').split(';') for line in a if 'KW  ' in line]
        keywords = [x.strip().lower() for list_ in KW_info for x in list_]
        keywords = [re.sub('.\{.*', '', key, flags=re.DOTALL) for key in keywords]
        keywords = [re.sub('.*\}.*', '', key, flags=re.DOTALL) for key in keywords] # remove weird codes
        keywords = list(filter(lambda x: x != "", keywords)) # check non empty strings
        if ECode:
            key = 'enzyme'
        elif any(key in class_list for key in keywords):
            for key in keywords:
                if key in class_list:
                    result = key; break
        else:
            key = 'notAnnotated' 
    data_uniprot.loc[len(data_uniprot.index)] = [uniprotid, SQ_info, key, ECode]



#### 
data_uniprot.to_pickle(os.path.join(OUT_PATH, 'data_uniprot_raw.pkl'))

frequencies = data_uniprot['Class'].value_counts()
condition = frequencies < 10
mask_obs = frequencies[condition].index
mask_dict = dict.fromkeys(mask_obs, 'notAnnotated')
data_uniprot['Class'] = data_uniprot['Class'].replace(mask_dict) 

#####################
data_uniprot.to_pickle(os.path.join(OUT_PATH, 'data_uniprot.pkl'))
# data_uniprot = pd.read_pickle(os.path.join(OUT_PATH, 'data_uniprot_raw.pkl'))

###########################################
# Quick check - results 

# import matplotlib.pyplot as plt
# import seaborn as sns
# df_full = data_uniprot

# df_full.UniprotID = df_full.UniprotID.astype(str) 
# df_crys = df_full[df_full.UniprotID.isin(list_proteins_struct)]

# sns.set_style("whitegrid")
# fig = plt.figure(figsize=(4,3),dpi=400)
# ax = fig.add_subplot(121)
# cts = df_full.Class.value_counts().to_frame()
# ax.pie(cts.Class)
# ax = fig.add_subplot(122)
# cts = df_crys.Class.value_counts().to_frame()
# ax.pie(cts.Class)
# #fig.legend(cts, fontsize='xx-small', loc="center right")
# plt.savefig('../Results/test_distr_test2.pdf')
# plt.clf()

# plt.title('crys')
# df_crys['Class'].value_counts().plot(kind='pie')
# plt.savefig('../Results/test_distr_test2.pdf')

# data.EC.fillna(value=np.nan, inplace=True)
# # select only first class
# data.EC  = data.EC.dropna().apply(lambda x: str(x)[0]).astype(int)
# #
# data = data.sort_values(by='Class').reset_index()
# data = data[data.Class == 'enzyme']
# mat_or = pd.DataFrame(matrix, columns=ordered, index=ordered)
# data_e = data.sort_values(by='EC')
# data_uniprot.EC.map(enz_dic)


