import os
import json
import pandas as pd
import pubchempy as pcp
import time
import numpy as np
import multiprocessing as mp

# ! CHANGE ABSOLUTE PATHS HERE

# this was used to optain the dic_drugnames_cid.json

output_path = '../Data/DrugBank'
data_path = '../../../Data/cross_side_information_DB/CTD'


# compound_id = str(pcp.get_compounds('10-decarbamoylmitomycin C', 'name')[0].cid) # class int -> str
def get_cid(x):
    try:
        comp_id = pcp.get_compounds(str(x), 'name')[0].cid
        time.sleep(1) # aumentar si sigue dando error
    except IndexError:
        comp_id = np.nan
        #print('ups, not in PubChem')
    return comp_id

###################################

h_chem = ['ChemicalName', 'ChemicalID' ,'CasRN', 'DiseaseName', 'DiseaseID',
        'DirectEvidence', 'InferenceGeneSymbol', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

chem_dis = pd.read_csv(os.path.join(data_path, 'CTD_chemicals_diseases.csv'), 
                            index_col=False, 
                            names=h_chem, 
                            comment='#',
                            usecols=['ChemicalName', 'DiseaseID'])

## remove NaN
chem_dis = chem_dis.dropna()
chem_dis.shape


# take all unique compounds, make a dictionary and then map
drugnames = chem_dis.ChemicalName.values.tolist()
drugnames = list(set(drugnames))

pool = mp.Pool(3) # multiprocss
cid_keys = pool.map(get_cid, drugnames)

dic_drugnames_cid = dict(zip(drugnames, cid_keys)) 

#save this dictionary in cross_side..._DB
file_name_json = '/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/dic_drugnames_cid.json'
with open(file_name_json, 'w', encoding='utf-8') as f:
    json.dump(dic_drugnames_cid, f, ensure_ascii=False, indent=4)
