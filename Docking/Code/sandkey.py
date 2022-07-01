
import os
import pandas as pd
import numpy as np
import helper_functions as hf
import logging


logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)



PKLS_PATH = '../Data/pkls'


# Annotations
drugs_annot = pd.read_pickle(os.path.join(PKLS_PATH, 'drugs_annot_full.pkl'))
drugs_annot.PubChemID = drugs_annot.PubChemID.astype(str)

# for Proteins
prot_annot = pd.read_pickle(os.path.join(PKLS_PATH, 'proteins_annot_crys.pkl'))


# Load all DTIS in a dicionary
dict_dfs  = {
             'DrugBank': hf.get_dtis_drugbank(),
             'Davis_et_al': hf.get_dtis_davis(),
             'BindingDB': hf.get_bindingdb_dtis(),
             'BIOSNAP': hf.get_biosnap_dtis()
             }

dict_dfs.update(hf.get_dict_dtis_yamanishi())


key = 'IC'
out_path = f'../Results/data_per_dataset/{key}'
logging.info(f'Working in {key}')
df = dict_dfs.get(key)

df['Protein_Family'] = df.UniprotID.map(dict(zip(prot_annot.id, prot_annot.Class )))
#df['Protein_Family'].value_counts()
df['Drug_Superclass'] = df.PubChemID.map(dict(zip(drugs_annot.PubChemID, drugs_annot.superclass)))
df = df.fillna('NotAnnotated')

#df.value_counts()
df.groupby(['Protein_Family']).size()
df.groupby(['Drug_Superclass']).size()


test_prot_ = pd.DataFrame(df.UniprotID.unique(), columns = ['UniprotID'])

test_prot_['l1'] =df.UniprotID.map(dict(zip(prot_annot.id, prot_annot.Class )))
test_prot_ = test_prot_.fillna('NotAnnotated')

test_prot_.groupby(['l1']).size()

df_chembl['uniprot'] = df_chembl.target_chembl_id.map(chembl2uni)

df['l1'] = df.UniprotID.map(dict(zip(df_chembl.uniprot, df_chembl.l1)))


test_prot = pd.DataFrame(df.UniprotID.unique(), columns = ['UniprotID'])

test_prot['l1'] = test_prot.UniprotID.map(dict(zip(df_chembl.uniprot, df_chembl.l1)))
df = test_prot.fillna('NotAnnotated')

test_prot.groupby(['l1']).size()

df_protein_annot.groupby('Class').size()