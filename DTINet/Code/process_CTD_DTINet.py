import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import multiprocessing as mp
import time
import json

# ! check absolute paths here

### IMPORTANT: This file uses all_proteins.txt which is created after running process_HPRD_DTINet.py

output_path = '../Data/DrugBank'
data_path = '../../../Data/cross_side_information_DB/CTD'


##############################################
####  Comparative Toxicogenomics Database ####

## Used in DTINet for:
##     - Disease nodes
##     - Drug-disease associations
##     - Protein-disease associations:

#########################################################################
############################ DISEASE NODES ##############################
# using: CTD_diseases.csv 

header_dis_dic = ['DiseaseName', 'DiseaseID', 'AltDiseaseIDs', 'Definition', 'ParentIDs', 
    'TreeNumbers', 'ParentTreeNumbers', 'Synonyms', 'SlimMappings']

dis_dict = pd.read_csv(os.path.join(data_path, 'CTD_diseases.csv'), 
                        index_col=False, 
                        names=header_dis_dic, 
                        comment='#', 
                        usecols=['DiseaseName', 'DiseaseID'])
dis_dict.DiseaseName = dis_dict.iloc[:, :1].applymap(lambda s: s.lower() if type(s) == str else s)

# In case we only want the number withouth MESH: or OMIM:
# dis_dict.DiseaseID = dis_dict.iloc[:, 1:2].applymap(lambda s: s[5:] if type(s) == str else s)


# export to files:::: diseases.txt is a list of diseases
np.savetxt(os.path.join(output_path,'disease.txt') , dis_dict.DiseaseName.values, newline='\n', fmt='%s')
# export to files:::: disease_dict_map.txt is a disctionary of diseases and identifiers
np.savetxt(os.path.join(output_path, 'disease_dict_map.txt'), dis_dict.DiseaseName.values + ',' + dis_dict.DiseaseID.values, newline='\n', fmt='%s')

dic_disease_to_disID = dict(list(zip(dis_dict.DiseaseName.values, dis_dict.DiseaseID.values)))
# or as json
file_name_json = 'dic_disease_to_disID.json'
with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
    json.dump(dic_disease_to_disID, f, ensure_ascii=False, indent=4)


###########################################################################################
##################### PROTEIN - DISEASE ASSOCIATIONS ######################################

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
gen_voc = gen_voc.explode('UniProtIDs') # as different entry


# now filter those that are not in the human protein set 
# extracted from HPRD (all_protein.txt) in output folder file

protein_file = os.path.join(output_path, 'all_protein.txt')
with open(protein_file, 'r') as fl:
    hum_prots = fl.read().splitlines()

## filter and keep only those rows that have a HUMAN UniProt ID !!
gen_voc = gen_voc.astype({"UniProtIDs": str})
gen_voc['comp'] = gen_voc['UniProtIDs'].isin(hum_prots)
gen_voc[gen_voc["comp"]==True].count()  # da 9061 uf
gene_and_uniprot = gen_voc[gen_voc["comp"]==True][['GeneSymbol','UniProtIDs']]

# create a dictinonary for mappig later to uniprot!
# genesymbol == keys ; uniprot ids == values
keys_gensymb = gene_and_uniprot['GeneSymbol'].values.tolist()
values_uniprot = gene_and_uniprot['UniProtIDs'].values.tolist()
dic_gen_to_protein_CTD = dict(list(zip(keys_gensymb,values_uniprot)))


#################
## genes diseases # use gene symbol as principal key (sq)
gen_dis_head = ['GeneSymbol', 'GeneID', 'DiseaseName', 'DiseaseID',
                'DirectEvidence', 'InferenceChemicalName', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

gen_dis = pd.read_csv(os.path.join(data_path, 'CTD_genes_diseases.csv'),
                        index_col=False, names=gen_dis_head, comment='#', usecols= ['GeneSymbol', 'DiseaseID'])

# gen_dis.shape
gen_dis.isnull().values.any()
gen_dis = gen_dis.dropna()
# gen_dis is a dataframe that relates each gene with a DiseaseID
# now use the dictionary created above to map
gen_dis['UniprotID'] = gen_dis['GeneSymbol'].map(dic_gen_to_protein_CTD)
gen_dis = gen_dis.dropna()
protein_disease = gen_dis[['UniprotID', 'DiseaseID']]
protein_disease = protein_disease.drop_duplicates()
protein_disease.to_csv(os.path.join(output_path, 'coordinates_protein_disease.tsv'), sep="\t")
# save like this and filter before creating the matix, 
# otherwise the matrix has ~5e7 rows 
######


##########################################################################
##################### DRUG - DISEASE ASSOCIATIONS ########################

h_chem = ['ChemicalName', 'ChemicalID' ,'CasRN', 'DiseaseName', 'DiseaseID',
        'DirectEvidence', 'InferenceGeneSymbol', 'InferenceScore', 'OmimIDs', 'PubMedIDs']

chem_dis = pd.read_csv(os.path.join(data_path, 'CTD_chemicals_diseases.csv'), 
                            index_col=False, 
                            names=h_chem, 
                            comment='#',
                            usecols=['ChemicalName', 'DiseaseID', 'DiseaseID'])

## remove NaN
chem_dis = chem_dis.dropna()
chem_dis.shape

# we need the DrugBankID, as is the identification that we use in DTINet
# Chemical ID is a MESH ID that appears in Pubchem, buy as synonym
# we cannot ask directly for a synonym to retrieve the CID
# then change from ChemicalName to PubchemCID using pubchempy

file_name_json = '/home/uveleiro/data/jfuente/DTI/Data/cross_side_information_DB/dic_drugnames_cid.json'

# read the dictionary (slow; just read but make the script above available) - get_dic_names_to_CID.py
with open(file_name_json, 'r') as f:
    dic_drugnames_cid = json.load(f)

# Then, we need to go back to DrugBank
# and make a relation between PubChem and DrugBank
################## esto esta repetido en SIDER ### se podrian juntar estos dos scripts
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
dic_cid_dbid[list(dic_cid_dbid.keys())[0]]
dic_cid_dbid.pop(list(dic_cid_dbid.keys())[0])


## apply two maps
chem_dis['PubChemID'] = chem_dis['ChemicalName'].map(dic_drugnames_cid)
chem_dis.shape
chem_dis = chem_dis.dropna()
chem_dis.shape # loosing 1158846 entries here
chem_dis['PubChemID'] = chem_dis['PubChemID'].astype(int).astype(str)

# esto funciona # usar lambda
# pp = np.unique(chem_dis['PubChemID'].values.tolist())
# tt =[x if x in pp else np.nan for x in dic_cid_dbid]

chem_dis['DrugBankID'] = chem_dis['PubChemID'].apply(lambda x: dic_cid_dbid[x] if x in dic_cid_dbid else np.nan)
chem_dis = chem_dis.dropna()

coordinates_drug_dis = chem_dis[['DrugBankID', 'DiseaseID']]
coordinates_drug_dis = coordinates_drug_dis.drop_duplicates()
coordinates_drug_dis.to_csv(os.path.join(output_path, 'coordinates_drug_dis.tsv'), index=False ,sep="\t")
# remove coordinates before creating matrix!
