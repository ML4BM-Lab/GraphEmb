import pandas as pd
import os
import seaborn as sns 
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib import pyplot as plt

#load dataframe
df = pd.read_csv(os.path.join('n2v_nn_results.tsv'),sep='\t',index_col = 0)
df['Embedding'] = df['Embedding'].apply(lambda x: int(x[1:]))
#drop train and validation results
df = df[['Network','Embedding','Architecture','Epochs','Loss Function','Test AUC','Test AUPRC']]

#separate auc and auprc
df_auc = df.drop('Test AUPRC',axis=1)
df_auprc = df.drop('Test AUC',axis=1)
df_dict = {'auc':df_auc, 'auprc':df_auprc}


"""
Violinplot
"""

nplots = len(set(df_auc['Network'].values))
networks = sorted(set(df_auc['Network'].values))
fig, axs = plt.subplots(nplots, 1, figsize=(15, 40))

for i in range(nplots):

    df_subplot = df_auc[df_auc['Network'] == networks[i]]
    df_subplot = df_subplot.sort_values(by='Embedding', ascending=False)

    #generate the violin plot
    sns.violinplot(x = "Embedding", y= "Test AUC", data = df_subplot, ax=axs[i], scale = 'width')
    axs[i].set_title(f'{networks[i]}')
    axs[i].set(ylim = (0,1.15))

#save as pdf
plt.savefig(os.path.join(f'n2v_nn_violinplot_splitted_Embedding_raw.pdf'))