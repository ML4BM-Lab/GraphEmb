import pandas as pd
import os
#!pip install PyTDC
from tdc.multi_pred import DTI

data = DTI(name = 'DAVIS')
data.label_distribution()
#data.convert_to_log(form = 'binding')
#data.convert_to_log(form = 'natural')
data.binarize(threshold = 5, order = 'descending')
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