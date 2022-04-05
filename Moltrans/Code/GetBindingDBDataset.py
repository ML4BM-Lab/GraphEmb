import pandas as pd
import os
#!pip install PyTDC
from tdc.multi_pred import DTI

data = DTI(name = 'BindingDB_Kd')
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