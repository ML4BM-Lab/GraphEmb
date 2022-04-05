# pip install PyTDC 

# duplicates already removed; data.harmonize_affinities(mode = 'max_affinity')
# https://tdcommons.ai/multi_pred_tasks/dti/#davis

from tdc.multi_pred import DTI


data = DTI(name = 'DAVIS')

data.binarize(threshold = 5, order = 'descending')
split = data.get_split()


len_sum = split['train'].shape[0] + split['test'].shape[0] + split['valid'].shape[0]
len_sum

test = split['test']
train = split['train']
valid = split['valid']
dataf = test.append(train)
dataf = dataf.append(valid)

dataf.shape[0] == len_sum

dataf.to_csv('Davis_5.tsv', sep="\t")

#####
from tdc.multi_pred import DTI

data = DTI(name = 'BindingDB_Kd')
data.harmonize_affinities(mode = 'max_affinity')