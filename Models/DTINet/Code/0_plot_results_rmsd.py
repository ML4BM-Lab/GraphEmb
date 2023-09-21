import matplotlib.pyplot as plt
import glob
import os
import re

import seaborn as sns
import pandas as pd

PATH = 'rmsd_v3_box_results'

files = glob.glob(os.path.join(PATH, 'BIOSNAP_*.out'))

#def get_data(file):

with open(file) as f:
    lines = f.readlines()


test_auroc_reg = '(?<=Test: AUROC=)[\d.]+'
test_aupr_reg = '(?<=AUPR=)[\d.]+(?=, FFold)'
ffold_auroc_reg = '(?<=FFold: AUROC=)[\d.]+'
ffold_aupr_reg = '(?<=AUPR=)[\d.]+$'

df = pd.DataFrame(columns = ['value', 'threshold', 'type'])

for file in files:
    threshold = int(re.search('(?<=BIOSNAP_)[\d]+', file).group())
    with open(file) as f:
        lines = f.readlines()
    for line in lines:
        df_line_t_auroc = {'value': float(re.search(test_auroc_reg, line).group()), 'threshold': threshold, 'type':'test_auroc'}
        df_line_t_aupr = {'value': float(re.search(test_aupr_reg, line).group()), 'threshold': threshold, 'type': 'test_aupr'}
        df_line_f_auroc =  {'value': float(re.search(ffold_auroc_reg, line).group()), 'threshold': threshold, 'type': 'bVal_auroc'}
        df_line_f_aupr =  {'value': float(re.search(ffold_aupr_reg, line).group()), 'threshold': threshold, 'type': 'bVal_aupr'}
        for d in [df_line_t_auroc, df_line_t_aupr, df_line_f_auroc, df_line_f_aupr]:
            df = df.append(d, ignore_index=True)


#list_types = ['val_auroc', 'val_aupr', 'test_auroc', 'test_aupr']

my_colors = ["#5c00b8","#C29CE7","#327ac3","#7CADDE"]
  
# add color array to set_palette
# function of seaborn
sns.set_palette( my_colors )

plt.clf()
plt.ylim(0.8,0.9)
plt.title('BIOSNAP v3 ')
plt.figure(figsize=(15, 8), dpi=300)
sns.pointplot(data=df, y='value', x='threshold', hue='type').set(title='BIOSNAP v3')
plt.xlabel("Threshold", size=20)
plt.ylabel("AUROC/AUPR", size=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.savefig('0_test_2.png', dpi=300)


#test = sns.load_dataset("penguins")