import pandas as pd
import os
import numpy as np
from scipy.stats import wilcoxon, mannwhitneyu
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns


## load reference
ref = pd.read_csv(os.path.join('RMSD_validation','final_fold_biosnap_annotated.tsv'), sep='\t', index_col=0)

## load --> final fold threshold rmsd 6
rmsd_files = sorted(os.listdir(os.path.join('RMSD_validation','FinalFold_threshold6_5_seeds_BIOSNAP')))
rmsd = {f'seed_{i}': pd.read_csv(os.path.join('RMSD_validation','FinalFold_threshold6_5_seeds_BIOSNAP',x,'predictions.csv'),index_col=0) for i,x in enumerate(rmsd_files)}

## load --> final fold threshold random
random_files = sorted(os.listdir(os.path.join('RMSD_validation','FinalFold_Random_5_seeds_BIOSNAP')))
random = {f'seed_{i}': pd.read_csv(os.path.join('RMSD_validation','FinalFold_Random_5_seeds_BIOSNAP', x, 'predictions.csv'),index_col=0) for i,x in enumerate(random_files)}


# pos_per_rmsd = []
# pos_per_random = []
# ## for each seed
# for i in range(5):
#     s = f'seed_{i}'

#     ## get only positives
#     only_positives = ref[ref['Label'] == 1].index

#     ## get pred and true
#     rmsd_score = rmsd[s].iloc[only_positives,0].values
#     rmsd_bool = rmsd[s].iloc[only_positives,1].values
#     pos_per_rmsd.append(rmsd_bool.sum()/rmsd_bool.shape[0])

#     random_score = random[s].iloc[only_positives,0].values
#     random_bool = random[s].iloc[only_positives,1].values
#     pos_per_random.append(random_bool.sum()/random_bool.shape[0])


# # perform wilcoxon test
# odds_ratio, p_value = wilcoxon(pos_per_rmsd, pos_per_random)

# # Create a DataFrame from the arrays
# data = pd.DataFrame({'RMSD': pos_per_rmsd, 'RANDOM': pos_per_random})
# # Plot boxplots using seaborn
# ax = sns.boxplot(data=data)
# ax.set_title(f'NS (p-val = {p_value})')
# # Add labels to the plot
# plt.xlabel('Conditions')
# plt.ylabel('Percentage')
# plt.savefig(os.path.join('RMSD_validation','rmsd_random_biosnap_comparison.png'))
# plt.clf()
# plt.cla()
# plt.close()


"""
CHECK SCORE

Gene	Protein_Name	UniProtKB	Drug ID	Drug_Name
P00533	Epidermal growth factor receptor	EGFR_HUMAN	DB08247	6-(cyclohexylmethoxy)-8-isopropyl-9h-purin-2-amine
P49841	Glycogen synthase kinase-3 beta	GSK3B_HUMAN	DB00748	carbinoxamine
"""

names = []
p1s = []
data = []

## get position of first
name = 'EGFR'
p1 = ref[ref['Gene'] == 'P00533']
p1 = p1[p1['Drug ID'] == 'DB08247'].index[0]

names.append(name)
p1s.append(p1)

name = 'GSKB'
p1 = ref[ref['Gene'] == 'P49841']
p1 = p1[p1['Drug ID'] == 'DB00748'].index[0]

names.append(name)
p1s.append(p1)

for name, p1 in zip(names,p1s):


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
    #wilcoxon(score_rmsd, score_random)
    odds_ratio, p_value = wilcoxon(score_rmsd, score_random)
    p_value = round(p_value,5)

    # Create a DataFrame from the arrays
    data.append(pd.DataFrame({'name':name, 'RMSD': score_rmsd, 'RANDOM': score_random}))


#build dataframe
df = pd.concat(data)
dfm = pd.melt(df, id_vars = ['name'])
dfm.columns = ['name','type','probability']

## set palette
custom_palette = ["#4f9ed5", "#5d5fa5"] 

# Set the custom palette
sns.set_palette(custom_palette)

# Plot boxplots using seaborn
ax = sns.boxplot(data=dfm, x = 'name', y = 'probability', hue='type')
ax.set_title(f'RMSD vs Random (5 seeds)')
# Add labels to the plot
plt.tight_layout()
plt.xlabel('Targets')
plt.ylabel('Prediction Probability')
plt.savefig(os.path.join('RMSD_validation',f'rmsd_random_biosnap_comparison_wilcoxon_final.pdf'))
plt.clf()
plt.cla()
plt.close()
