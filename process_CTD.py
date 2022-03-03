import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import multiprocessing as mp
import time
import json


####  Comparative Toxicogenomics Database ####

## Used in DTINet for:
##     - Disease nodes
##     - Drug-disease associations
##     - Protein-disease associations:

# saving outputs (test) in: /home/uveleiro/data/uveleiro/Test/output_test_CTD

# reading in
# os.chdir('/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/CTD')

# path to CTD folder

#########################################################################
############################ DISEASE NODES ##############################
# using: CTD_diseases.csv 

# also, this exist for retrieving info:  https://id.nlm.nih.gov/mesh/C538288.html

header_dis_dic = ['DiseaseName', 'DiseaseID', 'AltDiseaseIDs', 'Definition', 'ParentIDs', 
    'TreeNumbers', 'ParentTreeNumbers', 'Synonyms', 'SlimMappings']
dis_dict = pd.read_csv('/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/CTD/CTD_diseases.csv', 
                        index_col=False, 
                        names=header_dis_dic, 
                        comment='#', 
                        usecols=['DiseaseName', 'DiseaseID'])
dis_dict.DiseaseName = dis_dict.iloc[:, :1].applymap(lambda s: s.lower() if type(s) == str else s)
# In case we only want the number withouth MESH: or OMIM:
# dis_dict.DiseaseID = dis_dict.iloc[:, 1:2].applymap(lambda s: s[5:] if type(s) == str else s)


# export to files:::: diseases.txt is a list of diseases
np.savetxt('/home/uveleiro/data/uveleiro/Test/output_test_CTD/disease.txt', dis_dict.DiseaseName.values, newline='\n', fmt='%s')
# export to files:::: disease_dict_map.txt is a disctionary of diseases and identifiers
np.savetxt('/home/uveleiro/data/uveleiro/Test/output_test_CTD/disease_dict_map.txt', dis_dict.DiseaseName.values + ',' + dis_dict.DiseaseID.values, newline='\n', fmt='%s')



####################################################################
################# PROTEIN - DISEASE ASSOCIATIONS ###################

# gene csv # wget http://ctdbase.org/reports/CTD_genes_diseases.csv.gz

# read from CTD_genes.csv.

# gene vocabulary includes UniprotIDs
# but this is a mixed dataset, not only humans
# we can use the dictionaty created in HPRD, or a set from HPRD
# as the set contains only human protein UniprotKB IDs

header_gene_voc = ['GeneSymbol','GeneName','GeneID','AltGeneIDs','Synonyms','BioGRIDIDs','PharmGKBIDs','UniProtIDs']
gen_voc = pd.read_csv('/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/CTD/CTD_genes.csv',
                        index_col=False, names=header_gene_voc, comment='#', 
                        usecols=['GeneSymbol','UniProtIDs'])

gen_voc = gen_voc.dropna()
gen_voc.shape 

# los uniprots se pueden separar por '|' y ver si estan en la lista set proteis.txt
#gen_voc.UniProtIDs[35].split("|")
gen_voc.UniProtIDs = gen_voc.UniProtIDs.str.split("|")

# as different entry
gen_voc = gen_voc.explode('UniProtIDs')

# now filter those that are not in the human protein set 
# extracted from HPRD (proteins.txt)
protein_file = '/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/HPRD/proteins.txt'
with open(protein_file, 'r') as fl:
    hum_prots = fl.read().splitlines()

## filter and keep only those rows that have a HUMAN UniProt ID !!
gen_voc = gen_voc.astype({"UniProtIDs": str})
gen_voc['comp'] = gen_voc['UniProtIDs'].isin(hum_prots)
gen_voc[gen_voc["comp"]==True].count()  # da 9061 uf
gene_and_uniprot = gen_voc[gen_voc["comp"]==True][['GeneSymbol','UniProtIDs']]

# create a dictinonary for mappig later to uniprot!
# genesymbol == keys 
# uniprot ids == values
keys_gensymb = gene_and_uniprot['GeneSymbol'].values.tolist()
values_uniprot = gene_and_uniprot['UniProtIDs'].values.tolist()

# create the dictionary
dic_gen_to_protein_CTD = dict(list(zip(keys_gensymb,values_uniprot)))


#################
## genes diseases
# use gene symbol as principal key (sq)
gen_dis_head = ['GeneSymbol', 'GeneID', 'DiseaseName', 'DiseaseID',
                'DirectEvidence', 'InferenceChemicalName', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

gen_dis = pd.read_csv('/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/CTD/CTD_genes_diseases.csv',
                        index_col=False, names=gen_dis_head, comment='#', usecols= ['GeneSymbol', 'DiseaseID'])


# gen_dis.shape
gen_dis.isnull().values.any()
gen_dis = gen_dis.dropna()
# gen_dis.shape
# gen_dis is a dataframe that relates each gene with a DiseaseID
# now use the dictionary created above to map
gen_dis['UniprotID'] = gen_dis['GeneSymbol'].map(dic_gen_to_protein_CTD)
gen_dis = gen_dis.dropna()
protein_disease = gen_dis[['DiseaseID','UniprotID']]
protein_disease.to_csv('/home/uveleiro/data/uveleiro/Test/output_test_CTD/protein_disease.tsv', sep="\t")

##########################################################################
##########################################################################
##################### DRUG - DISEASE ASSOCIATIONS ########################

h_chem = ['ChemicalName', 'ChemicalID' ,'CasRN', 'DiseaseName', 'DiseaseID',
        'DirectEvidence', 'InferenceGeneSymbol', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

chem_dis = pd.read_csv('/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/CTD/CTD_chemicals_diseases.csv', 
                            index_col=False, 
                            names=h_chem, 
                            comment='#',
                            usecols=['ChemicalName', 'DiseaseID'])

## first clean the NaN but we only care about NaN in Chemical Name or in Disease ID
## indeed we can remove ChemicalID from the Dataset now
chem_dis.shape
chem_dis = chem_dis.dropna()
chem_dis.shape

# we need the DrugBankID, as is the identification that we use in DTINet
## Chemical ID is a MESH ID that appears in Pubchem, buy as synonym
## we cannot ask directly for a synonym to retrieve the CID
# then change from ChemicalName to PubchemCID using pubchempy
# -------> *******

# compound_id = str(pcp.get_compounds('10-decarbamoylmitomycin C', 'name')[0].cid) # class int -> str
def get_cid(x):
    try:
        comp_id = pcp.get_compounds(str(x), 'name')[0].cid
        time.sleep(1) # aumentar si sigue dando error
    except IndexError:
        comp_id = np.nan
        #print('ups, not in PubChem')
    return comp_id


#### ---------------------------------
# take all unique compounds
# make a dictionary and then map
drugnames = chem_dis.ChemicalName.values.tolist()
drugnames = list(set(drugnames))

# multiprocss
pool = mp.Pool(3)
cid_keys = pool.map(get_cid, drugnames)



dic_drugnames_cid = dict(zip(drugnames, cid_keys)) 

#save this dictionary in cross_side..._DB
file_name_json = '/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/dic_drugnames_cid.json'
with open(file_name_json, 'w', encoding='utf-8') as f:
    json.dump(dic_drugnames_cid, f, ensure_ascii=False, indent=4)

# read this json
with open(file_name_json, 'r') as f:
    data_dict = json.load(f)


##### dic_drugnames_cid[list(dic_drugnames_cid.keys())[0]]
##### data_dict['17531']
#####

# Then, we need to go back to DrugBank
# and make a relation between PubChem and DrugBank

import xml.etree.ElementTree as ET
tree = ET.parse('/home/margaret/data/jfuente/DTI/Data/cross_side_information_DB/DrugBank/Data/full_database.xml')
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




#### mapear series con dictionarios

#### ---------------------------------
## apply two maps
# first: dic_drugnames_cid
########## **********************************************
################################### ***** ojo porque habria que ver si saca 0s para comparar
### int('005') == 5 returns true, check before comparing
#### Instead of map(), apply and define a function to compare int() int()

chem_dis['PubChemID'] = chem_dis['ChemicalName'].map(dic_drugnames_cid)
# and later: dic_cid_dbid
chem_dis['DrugBankID'] = chem_dis['PubChemID'].map(dic_cid_dbid)

## remove nans
'''
# mejor no hacerlo de todo si nos interesan los pubchems
chem_dis.shape
chem_dis = chem_dis.dropna()
chem_dis.shape
'''

# -------> *******

## Change from Pubchem to CID


# --------> *******
# create a dictionary

# --------> *******


# for that we need to create a distionary to later map the Series
