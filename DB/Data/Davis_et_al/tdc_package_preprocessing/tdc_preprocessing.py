import numpy as np
import pandas as pd
import random
import os
#!pip install PyTDC
from tdc.multi_pred import DTI

data = DTI(name = 'DAVIS')
#data.harmonize_affinities(mode = 'max_affinity')
#data.harmonize_affinities(mode = 'mean')

#df = data.harmonize_affinities(mode = 'max_affinity') ONLY FOR BINDING DB
df = data.get_data() #FOR DAVIS
df = df.rename(columns = {'Drug':'SMILES', 
                        'Target':'Target Sequence'})

#Save
output_path = os.getcwd() + '/DAVIS_et_al.tsv'
df.to_csv(output_path,sep='\t')

#Load it
df = pd.read_csv(output_path,sep='\t',index_col=0)

drugs = df['SMILES'].unique()
genes = df['Target Sequence'].unique()

threshold = 30

df['Label'] = [1 if x < threshold else 0 for x in df['Y']]

print("Num Drugs:", len(drugs))
print("Num Genes:", len(genes))
print(df['Label'].value_counts())

#Save it as a csv files
df = pd.read_csv(output_path,sep='\t',index_col=0)

output_path_2 = os.getcwd() + '/DAVIS_et_al_w_labels.tsv'

df.to_csv(output_path_2,sep='\t')
