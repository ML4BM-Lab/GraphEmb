import pandas as pd
import os
#!pip install PyTDC
from tdc.multi_pred import DTI

data = DTI(name = 'DAVIS')
data.label_distribution()
data.binarize(threshold =30, order = 'descending')
data.label_distribution()

#harmonize affinities not supported because duplicities already removed
df = data.get_data()

df = df.rename(columns ={'Drug':'SMILES', 'Target':'Target Sequence'})

print(df.columns)

drugs = df['SMILES'].unique()
genes = df['Target Sequence'].unique()
print("Num Drugs:", len(drugs))
print("Num Genes:", len(genes))

#threshold = 30
#df['Label'] = [1 if x < threshold else 0 for x in df['Y']]

print(df['Y'].value_counts())

#Save it as a csv file
output_path = os.getcwd() + '/../Data/Davis_et_al/Davis.csv'
df.to_csv(output_path)


#Get data for splits
df_positives = df[df["Y"] == 1]
columns = df_positives[['Drug_ID','Target_ID', 'Y']]
df_splits = columns.copy()
output_path=os.getcwd() + '/../Data/Davis_et_al/Davis_pairs.txt'
df_splits.to_csv(output_path)


#Get drugs_smiles
df_drugs_smiles = df[['Drug_ID', 'SMILES']]
print(df_drugs_smiles)
df_drugs_smiles.drop_duplicates(inplace = True)
print(df_drugs_smiles)

output_path = os.getcwd() + '/../Data/Davis_et_al/drugs_smiles_Davis.txt'
df_drugs_smiles.to_csv(output_path, header=None, index = None, sep = ' ')

#Get gene target_sequences

df_genes_targetsequences = df[['Target_ID', 'Target Sequence']]
print(df_genes_targetsequences)
df_genes_targetsequences.drop_duplicates(inplace = True)
print(df_genes_targetsequences)

output_path = os.getcwd() + '/../Data/Davis_et_al/genes_targetsequences_Davis.txt'
df_genes_targetsequences.to_csv(output_path, header=None, index = None, sep = ' ')