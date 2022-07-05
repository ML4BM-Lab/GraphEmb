import pandas as pd
import os
#!pip install PyTDC
from tdc.multi_pred import DTI

data = DTI(name = 'BindingDB_Kd')
#data.convert_to_log(form = 'binding')
df = data.harmonize_affinities(mode = 'max_affinity')

#df = data.harmonize_affinities(mode = 'mean')
df = df.rename(columns ={'Drug':'SMILES', 'Target':'Target Sequence'})

drugs = df['SMILES'].unique()
genes = df['Target Sequence'].unique()

print("Num Drugs:", len(drugs))
print("Num Genes:", len(genes))

threshold = 30
df['Label'] = [1 if x < threshold else 0 for x in df['Y']]

print(df['Label'].value_counts())

#Save it as a csv files
output_path = os.getcwd() + '/../Data/BindingDB/BindingDB.csv'
df.to_csv(output_path)

#Get data for splits
df_positives = df[df["Label"] == 1]
print(df_positives)
columns = df_positives[['Drug_ID','Target_ID', 'Label']]
df_splits = columns.copy()
output_path=os.getcwd() + '/../Data/BindingDB/BindingDB_pairs.txt'
df_splits.to_csv(output_path)

#Get drugs_smiles
df_drugs_smiles = df_positives[['Drug_ID', 'SMILES']]
print(df_drugs_smiles)
df_drugs_smiles.drop_duplicates(inplace = True)
print(df_drugs_smiles)

output_path = os.getcwd() + '/../Data/BindingDB/drugs_smiles_BindingDB.txt'
df_drugs_smiles.to_csv(output_path, header=None, index = None, sep = ' ')

#Get gene target_sequences

df_genes_targetsequences = df_positives[['Target_ID', 'Target Sequence']]
print(df_genes_targetsequences)
df_genes_targetsequences.drop_duplicates(inplace = True)
print(df_genes_targetsequences)

output_path = os.getcwd() + '/../Data/BindingDB/genes_targetsequences_BindingDB.txt'
df_genes_targetsequences.to_csv(output_path, header=None, index = None, sep = ' ')