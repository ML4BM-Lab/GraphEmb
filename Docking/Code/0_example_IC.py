import pandas as pd
import os

PKLS_PATH = '../Data/pkls'

prot_rmsd = pd.read_pickle(os.path.join(PKLS_PATH, 'RMSD_full_matrix.pkl')) #
prot_annot = pd.read_pickle(os.path.join(PKLS_PATH, 'proteins_annot_crys.pkl'))

df_chembl = None # loaded from get_protein_annot_chembl
# chembl2uni dict
prots = None # list from IC

df_chembl['UniprotID'] = df_chembl.target_chembl_id.map(chembl2uni) 


prots

df = pd.DataFrame()
df['UniprotID'] = prots


df['from_uniprot'] = df['UniprotID'].map(dict(zip(prot_annot.id, prot_annot.Class)))
df['from_chembl'] = df['UniprotID'].map(dict(zip(df_chembl.UniprotID, df_chembl.l1)))

df_test = df[~df.from_uniprot.isna()]
df_test = df_test[~df_test.from_chembl.isna()]


df_test.from_chembl = df_test.from_chembl.apply(lambda x: x.lower())

df_test_ = df_test[df_test.from_uniprot != df_test.from_chembl]
df_test_.to_excel('.test_IC_Data_meeting.xlsx')
