import os
import pandas as pd


#read df
# bindingdb = pd.read_csv(os.path.join('.','BindingDB_max_affinity.tsv'),sep='\t')

# bindingdb_df = bindingdb.loc[:,['Drug_ID','Target_ID','Y']]

# threshold = 30
# bindingdb_df['Label'] = [1 if x < threshold else 0 for x in bindingdb_df['Y']]

# bindingdb_df.to_csv(os.path.join('.','BindingDB_max_affinity_w_labels.tsv'),sep='\t')

#bindingdb_df.to_csv(os.path.join('.','BindingDB_max_affinity_w_labels.tsv'),sep='\t')

bindingdb_df = pd.read_csv(os.path.join('.','BindingDB_max_affinity_w_labels.tsv'),sep='\t',index_col=0)

#get only 2 cols
bindingdb2line = bindingdb_df[bindingdb_df.Label == 1].loc[:,['Drug_ID','Target_ID']]


bindingdb2line.to_csv(os.path.join('.','BindingDB_max_affinity_w_labels_2line.tsv'),sep='\t',header=False,index=False)