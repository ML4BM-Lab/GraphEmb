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

with open(os.path.join('figure2','all_results_5fold.pickle'), 'rb') as handle:
    all_results = pickle.load(handle)

#build the proper df
all_results_df = pd.DataFrame(all_results)
all_results_df = pd.concat([all_results_df,all_results_df.iloc[:,2].apply(pd.Series)], axis=1)
all_results_df.columns = range(all_results_df.shape[1])
all_results_df = all_results_df.drop([2,6,7],axis=1)
all_results_df.columns= ['train_data','test_data','seed','train_auc','test_auc']

train_auc = pd.DataFrame(all_results_df.groupby(['train_data','test_data'])['train_auc'].mean()).reset_index()
test_auc = pd.DataFrame(all_results_df.groupby(['train_data','test_data'])['test_auc'].mean()).reset_index()

datasets = ['DrugBank','BIOSNAP','BindingDB','Davis_et_al','Yamanishi/E', 'Yamanishi/IC', 'Yamanishi/GPCR', 'Yamanishi/NR']
train_auc = train_auc.pivot(index='train_data',columns='test_data',values='train_auc').loc[datasets,datasets]
test_auc = test_auc.pivot(index='train_data',columns='test_data',values='test_auc').loc[datasets,datasets]


for var in ['train_auc','test_auc']:

    #divide in four
    heatmap_matrix = eval(var)

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
    data = heatmap_matrix.values.flatten()
    heatmap_matrix_norm = (heatmap_matrix - data.min()) / (data.max() - data.min())


    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(heatmap_matrix_norm, cmap=cmap, square=True, linewidths = 0.2, linecolor = 'gray', annot = heatmap_matrix, annot_kws={"size": 12})
    plt.savefig(os.path.join(f'heatmap_{var}.pdf'))