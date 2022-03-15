import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import json
import xml.etree.ElementTree as ET



def db_m(x):
    if x['Pubchem_map_flat'] == x['Pubchem_map_stereo']:
        return x['Pubchem_map_flat']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and ( x['Pubchem_map_stereo']) != np.nan ) :
        return x['Pubchem_map_stereo']
    elif (x['Pubchem_map_flat'] !=x['Pubchem_map_stereo'] and (np.isnan(x['Pubchem_map_flat']) == False)):
        return x['Pubchem_map_flat']
    else:
        return np.nan


###########
###########


tree = ET.parse(data_path_DrugBank)
root = tree.getroot()

dbids = []
pubchemids = []
for drug_entry in tqdm(root):
    drugbank_ID = drug_entry.find('{http://www.drugbank.ca}drugbank-id').text
    pubchem_id = np.nan # to not repeat id if it does not appear later
    for props in drug_entry.findall('.//{http://www.drugbank.ca}external-identifier'):
        for prop in props:
            if(prop.text == 'PubChem Compound'): 
                pubchem_id = props[1].text
            break # romper una vez que encuentre el pubchem 
    dbids.append(drugbank_ID)
    pubchemids.append(pubchem_id)

dic_cid_dbid = dict((zip(pubchemids, dbids)))
# first element is nan, delete
dic_cid_dbid[list(dic_cid_dbid.keys())[0]]
dic_cid_dbid.pop(list(dic_cid_dbid.keys())[0])

'''
# pubchemids as str here
with open(os.path.join(output_path ,'dic_cid_dbid.json'), 'w', encoding='utf-8') as f:
    json.dump(dic_cid_dbid, f, ensure_ascii=False, indent=4)
'''

##########
##########

output_path = '../Data/DrugBank'
data_path = '../../../Data/cross_side_information_DB/SIDER'
data_path_DrugBank = '/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml'

file = 'meddra_all_se.tsv'
df_drug_se = pd.read_csv(os.path.join(data_path, file), sep='\t',
            names = ['STITCH_flat', 'STITCH_stereo', 'x1', 'x2', 'x3', 'se'],
            usecols = ['se', 'STITCH_flat', 'STITCH_stereo'])

# lowerase
df_drug_se.se = df_drug_se.iloc[:, 2:3].applymap(lambda s: s.lower() if type(s) == str else s)
se_list = df_drug_se.se.unique().tolist()

# save se.txt
np.savetxt(os.path.join(output_path ,'se.txt'), se_list , newline='\n', fmt='%s')



#################

# key_list = list(dic_cid_dbid.keys())
# key_list_int = list(map(int, key_list)) #str

## take both to get the max posible number of DBIDs
df_drug_se['Pubchem_flat'] = df_drug_se['STITCH_flat'].apply(lambda x: x[4:]).astype(int)
df_drug_se['Pubchem_stereo'] = df_drug_se['STITCH_stereo'].apply(lambda x: x[4:]).astype(int)

df_drug_se['Pubchem_flat'] = df_drug_se['Pubchem_flat'].astype(str)
df_drug_se['Pubchem_stereo'] = df_drug_se['Pubchem_stereo'].astype(str)

df_drug_se['Pubchem_map_flat'] = df_drug_se['Pubchem_flat'].map(dic_cid_dbid)
df_drug_se['Pubchem_map_stereo'] = df_drug_se['Pubchem_stereo'].map(dic_cid_dbid)


# clean the dataset and we hace this!
#data_drug_to_se = df_drug_se[['se', 'se_ID', 'Pubchem_map_flat', 'Pubchem_map_stereo']]
data_drug_to_se = df_drug_se[['se', 'Pubchem_map_flat', 'Pubchem_map_stereo']].copy(deep=True)

# below the column to fill

data_drug_to_se['DrugBank_ID'] = data_drug_to_se.apply(db_m, axis=1)
data_drug_to_se = data_drug_to_se[['se', 'DrugBank_ID']]
data_drug_to_se = data_drug_to_se.dropna()
data_drug_to_se

coordinates_drug_se = data_drug_to_se[['se','DrugBank_ID']]
coordinates_drug_se.shape
coordinates_drug_se = coordinates_drug_se.drop_duplicates()
coordinates_drug_se.shape
coordinates_drug_se.to_csv(os.path.join(output_path, 'coordinates_drug_se.tsv'), header=True, index=False ,sep="\t")




df_luo = pd.read_csv(os.path.join(output_path, 'tmp/drug.txt'), header=None)

len(set(np.unique(df_luo.values)).intersection(set(coordinates_drug_se.DrugBank_ID.unique())))