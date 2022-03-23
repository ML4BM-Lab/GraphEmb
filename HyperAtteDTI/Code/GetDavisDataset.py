import pandas as pd
import numpy as np
import os
from pubchempy import Compound
from collections import Counter
import requests

from tdc.multi_pred import DTI
data = DTI(name = 'DAVIS')

#data.binarize(threshold = 30, order = 'descending')

df = data.get_data()
#print(df.columns)
df = df.rename(columns ={'Drug':'SMILES', 'Target':'Target Sequence'})

drugs = df['SMILES'].unique()
genes = df['Target Sequence'].unique()

df['Label'] = [1 if x < 30 else 0 for x in df['Y']]
#split = data.get_split()

print("Num Drugs:", len(drugs))
print("Num Genes:", len(genes))

print(df['Label'].value_counts())

print(df.columns)

cols = ['Drug_ID', 'Target_ID', 'SMILES', 'Target Sequence', 'Label']
df = df[cols]

output_path = os.getcwd() + '/../Data/Davis_et_al/Davis.txt'
df.to_csv(output_path, header=None, index = None, sep = ' ')
#print(split)