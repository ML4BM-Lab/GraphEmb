from cmath import nan
import os
import numpy as np
import pandas as pd


########### BUILD A DICTIONARY ##################
## HPRD needs a dictorary to convert from genesymbol to uniprot
# so we can know which proteins are we using in our model

# first: generate  genes_data.dat  --> execute -> get_genes_data.sh

data = pd.read_table('../../Data/cross_side_information_DB/HPRD/genes_data.dat', names=['Uniprot_ID', 'type', 'GeneSymbol'])
translator = data.drop_duplicates(subset='GeneSymbol', keep="first") # GeneName appears always before GeneCards, so its safe
# translator[translator['GeneSymbol']=="YWHAB"] # check appears once
# get value list for dictionary
val_protein_id = translator['Uniprot_ID'].values.tolist()
key_gen_name = translator['GeneSymbol'].values.tolist()
# create the dictionary
dic_gen_to_protein = dict(list(zip(key_gen_name,val_protein_id)))
# dic_gen_to_protein.get('YWHAB', 0) #check first element



################# GET HPDR DATA ###############
######### PROTEIN - PROTEIN INTERACTIONS ######
hpdr_file = '../../Data/cross_side_information_DB/HPRD/HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'

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


hprd_processed.shape
hprd_processed.head()

#dic_gen_to_protein.get("ETS1",None)


path = '../../Data/cross_side_information_DB/HPRD/'
hprd_processed.to_csv(os.path.join(path,'hprd_processed_with_uniprot.tsv'), sep="\t")


###########################################################################
########################### PROTEIN NODES #################################

s1 = hprd_processed['Prot1 UniprotKB'].values.tolist()
s2 = hprd_processed['Prot2 UniprotKB'].values.tolist()
s = list(set(s1 + s2))
# list of proteins in HPRD (implies Human)
np.savetxt(os.path.join(path, 'proteins.txt'), s, newline='\n', fmt='%s')
