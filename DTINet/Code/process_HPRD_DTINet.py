import os
from re import M
import numpy as np
import pandas as pd
import json
import requests
from tqdm import tqdm

### from HPRD we need
#       - protein.txt 
#       - protein_dict_map.txt
#       - mat_protein_protein.txt
#    ****  - Similarity_Matrix_Proteins.txt (Smith-Waterman score) -> (* script get_SW_score find) **** 

output_path = '../Data/DrugBank'
data_path = '../../../Data/cross_side_information_DB/HPRD'

####################################################################
############### BUILD A DICTIONARY GEN-UNIPROTID ###################

# for HPRD we need a dictorary to convert from genesymbol to uniprot

# first: generate  genes_data.dat  --> execute -> get_genes_data.sh
os.system(f'grep "Gene_Name\|GeneCards" {data_path}/HUMAN_9606_idmapping.dat > {output_path}/genes_data.dat')

''' CHANGED!!!! for this : 
with open(os.path.join(data_path,'HUMAN_9606_idmapping.dat' ), 'r') as infl:
    all_lines = infl.read().splitlines()

all_lines = list(filter(lambda line: ('Gene_Name'  in line) or ('GeneCards' in line), all_lines))
all_lines_spt = [i.split('\t') for i in all_lines] 
test_dt = pd.DataFrame(all_lines_spt, columns=['Uniprot_ID', 'type', 'GeneSymbol'])
'''
# test_dt.equals(data) IS TRUE ==> replace !!

# get the data for processing
data = pd.read_table(os.path.join(output_path, 'genes_data.dat'), names=['Uniprot_ID', 'type', 'GeneSymbol'])
translator = data.drop_duplicates(subset='GeneSymbol', keep="first") # GeneName appears always before GeneCards, so its safe
# translator[translator['GeneSymbol']=="YWHAB"] # check appears once
# get value list for dictionary
val_protein_id = translator['Uniprot_ID'].values.tolist()
key_gen_name = translator['GeneSymbol'].values.tolist()
# create the dictionary
dic_gen_to_protein = dict(list(zip(key_gen_name,val_protein_id)))
# dic_gen_to_protein.get('YWHAB', 0) #check first element

#### The dictionary that they have is:
dic_protein_to_gen = dict(list(zip(val_protein_id,key_gen_name)))

# Write dictionary (txt file) (drugbank_ID and name)
with open(os.path.join(output_path,'protein_dic_map.txt'), 'w') as f:
	for i in range(len(val_protein_id)):
		_ = f.write("%s:%s\n" % (val_protein_id[i],key_gen_name[i]))

# --------- uncomment to save as json -------------
#file_name_json = 'dic_protein_to_gen.json'
#with open(os.path.join(output_path,file_name_json), 'w', encoding='utf-8') as f:
#    json.dump(dic_protein_to_gen, f, ensure_ascii=False, indent=4)


####################################################################
################## PROTEIN - PROTEIN INTERACTIONS ##################

hpdr_file = os.path.join(data_path, 'HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt')

hpdr_head = ['Prot1 GeneSymbol', 'Prot1 HPRDid', 'Prot1 RefSeqid',
    'Prot2 GeneSymbol', 'Prot2 HPRDid', 'Prot2 RefSeqid',
    'Experiment type', 'Pubmed id']
hprd = pd.read_table(hpdr_file, names=hpdr_head, usecols=['Prot1 GeneSymbol', 'Prot1 HPRDid','Prot2 GeneSymbol', 'Prot2 HPRDid'])

hprd['Prot1 UniprotKB'] = hprd['Prot1 GeneSymbol'].map(dic_gen_to_protein)
hprd['Prot2 UniprotKB'] = hprd['Prot2 GeneSymbol'].map(dic_gen_to_protein)
# drop HPDR (dont actually need tat)
hprd_processed =  hprd.drop(columns=['Prot1 HPRDid','Prot2 HPRDid'])
# first drop rows that contain no Gene Symbol ('-')
hprd_processed = hprd_processed.drop(index=hprd_processed[hprd_processed['Prot1 GeneSymbol']=='-'].index, axis=1)
hprd_processed = hprd_processed.drop(index=hprd_processed[hprd_processed['Prot2 GeneSymbol']=='-'].index, axis=1)
# drop nan
hprd_processed = hprd_processed.dropna()
hprd_processed[hprd_processed['Prot1 UniprotKB'].isna()]
hprd_processed[hprd_processed['Prot2 UniprotKB'].isna()]

# can this can be saved here as
# hprd_processed.to_csv(os.path.join(output_path,'hprd_processed_with_uniprot.tsv'), sep="\t")

# Create a squared matrix of protein-protein interaction
df_PPI = hprd_processed.drop(['Prot1 GeneSymbol', 'Prot2 GeneSymbol'], axis=1)
df_PPI.columns = ['Prot1', 'Prot2']
df_PPI = df_PPI.drop_duplicates()
# save coordinate files
df_PPI.to_csv(os.path.join(output_path, 'coordinates_PPI.tsv'), index=False, header=False, sep='\t') 

### matrix
'''
# last step; now get only coordinates! 
matrix_prot_prot_ = pd.get_dummies(df_PPI.set_index('Prot1')['Prot2']).max(level=0)#.reset_index()# no esta ordenado esto
matrix_prot_prot_.columns
# their matrix is also squared
matrix_prot_prot = pd.DataFrame(matrix_prot_prot_, columns= matrix_prot_prot_.columns, index= matrix_prot_prot_.columns)
matrix_prot_prot = matrix_prot_prot.fillna(0)
matrix_prot_prot.head(2)
# save
matrix_prot_prot.to_csv(os.path.join(output_path, 'mat_protein_protein.tsv'), index=True, header=True, sep='\t') 
'''

###########################################################################
########################### PROTEIN NODES #################################

s1 = hprd_processed['Prot1 UniprotKB'].values.tolist()
s2 = hprd_processed['Prot2 UniprotKB'].values.tolist()
protein_list = list(set(s1 + s2))
# list of proteins in HPRD (implies Human)
np.savetxt(os.path.join(output_path, 'all_protein.txt'), protein_list, newline='\n', fmt='%s')