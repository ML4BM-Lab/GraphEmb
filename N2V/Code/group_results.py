import os
os.environ["OMP_NUM_THREADS"] = "16"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import json
import pandas as pd
import itertools
import numpy as np
from collections import Counter
import importlib
import pickle
import sklearn.metrics
from tqdm import tqdm
import sys
from collections import defaultdict as dd


## Detect the different networks
n2v_nn_results = os.listdir(os.path.join('n2v_nn_results'))
#sorted set n2v_nn results
n2v_nn_results_ss = sorted(set(['_'.join(x.split('_')[0:2]) if x.split('_')[0] == 'Yamanishi' else x.split('_')[0] for x in n2v_nn_results]))

## Define a network dict
network_dd = dd(list)

for n2v_nn_result in n2v_nn_results:

    #get which network is this pickle coming from
    network = n2v_nn_results_ss[np.where([x in n2v_nn_result for x in n2v_nn_results_ss])[0][0]]

    #append to dict
    network_dd[network].append(n2v_nn_result)


#define pickles path
n2v_nn_results_path = os.path.join('n2v_nn_results')


full_pickle_df = []
for network in network_dd:

    #read all pickles
    pickle_l = [[network] + pd.read_pickle(os.path.join(n2v_nn_results_path, pickle_n2v_nn)) for pickle_n2v_nn in network_dd[network]]
    
    #convert to df
    pickle_df = pd.DataFrame(pickle_l, columns=['Network', 'Train AUC', 'Val AUC', 'Test AUC','Train AUPRC','Val AUPRC', 'Test AUPRC'], index=network_dd[network])
    
    #concat
    full_pickle_df.append(pickle_df)
    
    #sort and print first 5
    print(pickle_df.sort_values(by=['Test AUPRC'], ascending=False).head(5)) 


def get_params(_str):

    #trim it
    tstr = _str[:-7]
    #split it
    ststr = tstr.split('_')
    loss_f = ststr[-1]
    arc = ' '.join(ststr[-3:-1])
    epochs = ststr[-4]
    emb = ststr[-9]

    return emb, arc, epochs, loss_f

#concat and add params
full_df = pd.concat(full_pickle_df)
params_df = pd.DataFrame(list(map(get_params, full_df.index)),
                        columns = ['Embedding','Architecture','Epochs', 'Loss Function'], 
                        index = full_df.index)

#save
full_df_params = pd.concat([params_df,full_df],axis=1)
full_df_params_order = full_df_params[['Network', 'Embedding', 'Architecture', 'Epochs',
                                       'Loss Function', 'Train AUC', 'Val AUC', 'Test AUC', 
                                       'Train AUPRC', 'Val AUPRC', 'Test AUPRC']]
full_df_params_order.to_csv(os.path.join('panel','figure1','n2v_nn_results.tsv'),sep='\t')