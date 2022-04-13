import pandas as pd
import os
#!pip install PyTDC
from tdc.multi_pred import DTI

data = DTI(name = 'DAVIS')

#harmonize affinities not supported because duplicities already removed
df = data.get_data()

df = df.rename(columns ={'Drug':'SMILES', 'Target':'Target Sequence'})

drugs = df['SMILES'].unique()
genes = df['Target Sequence'].unique()
print("Num Drugs:", len(drugs))
print("Num Genes:", len(genes))

threshold = 30
df['Label'] = [1 if x < threshold else 0 for x in df['Y']]

print(df['Label'].value_counts())

#Save it as a csv file
output_path = os.getcwd() + '/../Data/Davis_et_al/Davis.csv'
df.to_csv(output_path)