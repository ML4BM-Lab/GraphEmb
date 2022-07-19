import helper_functions as hf
import pandas as pd

# Load all DTIS in a dicionary
dict_dfs  = {
             'DrugBank': hf.get_dtis_drugbank(),
             'Davis_et_al': hf.get_dtis_davis(),
             'BindingDB': hf.get_bindingdb_dtis(),
             'BIOSNAP': hf.get_biosnap_dtis()
             }

dict_dfs.update(hf.get_dict_dtis_yamanishi())


key = 'DrugBank'
out_path = f'../Results/data_per_dataset/{key}'
df = dict_dfs.get(key)
df.groupby(by='PubChemID').count()

import matplotlib.pyplot as plt
plt.clf()
df.value_counts('PubChemID').plot.hist(bins = 100) 
plt.savefig('test_count.png')