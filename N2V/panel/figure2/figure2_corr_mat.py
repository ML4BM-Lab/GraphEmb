import pickle
import pandas as pd
import os
import seaborn as sns

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt

import numpy as np

def line_to_mat(lines):

    #define mat
    m = lines.shape[0]
    mroot = int(m**(1/2))

    #define mat
    mat = np.zeros(shape = (mroot,mroot))

    for i in range(mroot):
        for j in range(mroot):
            mat[i,j] = lines.iloc[i*mroot+j,2]

    df = pd.DataFrame(mat)
    df.columns = lines['test'].values[0:mroot]
    df.index = lines['test'].values[0:mroot]
    
    return df

# with open(os.path.join('figure2','only_diagonal_all_results.pickle'), 'rb') as handle:
#     all_results = pickle.load(handle)

with open(os.path.join('figure2','all_results.pickle'), 'rb') as handle:
    all_results = pickle.load(handle)

#build the proper df
all_results_df = pd.DataFrame(all_results)
multiple_columns = all_results_df.iloc[:,2].apply(pd.Series)

#concat
all_results_df_concat = pd.concat([all_results_df.iloc[:,[0,1]],multiple_columns],axis=1)
all_results_df_concat.columns = ['train','test','train_auc','test_auc','train_auprc','test_auprc']

for var in ['train_auc','test_auc']:

    #divide in four
    corr = all_results_df_concat.loc[:,['train','test',var]]
    corr.drop_duplicates(inplace=True)

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(20, 20))

    # Generate a custom diverging colormap
    #cmap = sns.diverging_palette(230, 20, as_cmap=True)
    if 'train' in var:
        #cmap = sns.color_palette("YlOrBr", as_cmap=True)
        cmap = sns.color_palette("icefire", 20,  as_cmap=True).reversed()
    else:
        #cmap = sns.color_palette("icefire", 20,  as_cmap=True).reversed()
        cmap = sns.color_palette("viridis", 20, as_cmap=True)
        #cmap = sns.color_palette("icefire",20,  as_cmap=True).reversed()

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(line_to_mat(corr), cmap=cmap, center=0, square=True, linewidths = 0.2)
    plt.savefig(os.path.join(f'heatmap_{var}.pdf'))