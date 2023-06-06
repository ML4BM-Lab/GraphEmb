import pandas as pd
import os
import numpy as np
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns


## load reference
ref = pd.read_csv(os.path.join('RMSD','final_fold_biosnap_annotated.tsv'), sep='\t', index_col=0)

## load --> final fold threshold rmsd 6
rmsd_files = sorted(os.listdir(os.path.join('RMSD','FinalFold_threshold6_5_seeds_BIOSNAP')))
rmsd = {f'seed_{i}': pd.read_csv(os.path.join('RMSD','FinalFold_threshold6_5_seeds_BIOSNAP',x,'predictions.csv'),index_col=0) for i,x in enumerate(rmsd_files)}

## load --> final fold threshold random
random_files = sorted(os.listdir(os.path.join('RMSD','FinalFold_Random_5_seeds_BIOSNAP')))
random = {f'seed_{i}': pd.read_csv(os.path.join('RMSD','FinalFold_Random_5_seeds_BIOSNAP', x, 'predictions.csv'),index_col=0) for i,x in enumerate(random_files)}


pos_per_rmsd = []
pos_per_random = []
## for each seed
for i in range(5):
    s = f'seed_{i}'

    ## get only positives
    only_positives = ref[ref['Label'] == 1].index

    ## get pred and true
    rmsd_score = rmsd[s].iloc[only_positives,0].values
    rmsd_bool = rmsd[s].iloc[only_positives,1].values
    pos_per_rmsd.append(rmsd_bool.sum()/rmsd_bool.shape[0])

    random_score = random[s].iloc[only_positives,0].values
    random_bool = random[s].iloc[only_positives,1].values
    pos_per_random.append(random_bool.sum()/random_bool.shape[0])


# perform wilcoxon test
odds_ratio, p_value = wilcoxon(pos_per_rmsd, pos_per_random)

# Create a DataFrame from the arrays
data = pd.DataFrame({'RMSD': pos_per_rmsd, 'RANDOM': pos_per_random})
# Plot boxplots using seaborn
ax = sns.boxplot(data=data)
ax.set_title(f'NS (p-val = {p_value})')
# Add labels to the plot
plt.xlabel('Conditions')
plt.ylabel('Percentage')
plt.savefig(os.path.join('RMSD','rmsd_random_biosnap_comparison.png'))
plt.clf()
plt.cla()
plt.close()


"""
CHECK SCORE

Gene	Protein_Name	UniProtKB	Drug ID	Drug_Name
P00533	Epidermal growth factor receptor	EGFR_HUMAN	DB08247	6-(cyclohexylmethoxy)-8-isopropyl-9h-purin-2-amine
P49841	Glycogen synthase kinase-3 beta	GSK3B_HUMAN	DB00748	carbinoxamine
"""

## get position of first
name = 'EGFR'
p1 = ref[ref['Gene'] == 'P00533']
p1 = p1[p1['Drug ID'] == 'DB08247'].index[0]

# name = 'GSKB'
# p1 = ref[ref['Gene'] == 'P49841']
# p1 = p1[p1['Drug ID'] == 'DB00748'].index[0]

score_rmsd = []
bool_rmsd = []
score_random = []
bool_random = []

for i in range(5):
    s = f'seed_{i}'

    ## get pred and true
    rmsd_score = rmsd[s].iloc[p1,0]
    rmsd_bool = rmsd[s].iloc[p1,1]
    score_rmsd.append(rmsd_score)
    bool_rmsd.append(rmsd_bool)

    random_score = random[s].iloc[p1,0]
    random_bool = random[s].iloc[p1,1]
    score_random.append(random_score)
    bool_random.append(random_bool)


# perform wilcoxon test
odds_ratio, p_value = wilcoxon(score_rmsd, score_random)

# Create a DataFrame from the arrays
data = pd.DataFrame({'RMSD': score_rmsd, 'RANDOM': score_random})
# Plot boxplots using seaborn
ax = sns.boxplot(data=data)
ax.set_title(f'{name}: NS (p-val = {p_value})')
# Add labels to the plot
plt.xlabel('Conditions')
plt.ylabel('Prediction Probability')
plt.savefig(os.path.join('RMSD',f'rmsd_random_biosnap_comparison_{name}.png'))
plt.clf()
plt.cla()
plt.close()
