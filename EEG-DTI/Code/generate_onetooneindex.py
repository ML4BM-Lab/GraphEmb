import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from tqdm import tqdm
import logging
import os
### we need 01
# => this is, protein drug

# load mat_protein_drug.txt

'''
logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)
'''

# test this now!
path = '../Data/Davis_et_al/'
folder = 'oneTooneIndex'
# still need to get this one in the previous step
# df = pd.read_csv('mat_drug_protein.txt', sep=" ", header=None)
# prote_drug = df.T
# prote_drug.to_csv('mat_protein_drug.txt', index=False, header=False, sep=" ")

admat = pd.read_csv(os.path.join(path, 'sevenNets/mat_protein_drug.txt'), sep=' ', header=None)
admat.shape

# select edges
rows = admat.index.tolist()
columns = admat.columns.tolist()
edges_all_ = []
edges_all_false_ = []
# iterate
for row in tqdm(rows):
    #print(row, end='\r')
    for column in columns:
        element = admat.loc[row][column]
        if element == 1:
            onetoone = (rows.index(row), columns.index(column))
            #print(onetoone)
            edges_all_.append(onetoone)
        elif element == 0:
            onetoone = (rows.index(row), columns.index(column))
            edges_all_false_.append(onetoone)
        else:
            print('edge nor 0 or 1!!')


edges_all = pd.DataFrame(edges_all_)
logging.info(f'# edges postives: {edges_all.shape[0]} ')
edges_all_false = pd.DataFrame(edges_all_false_)
logging.info(f'# edges false: {edges_all_false.shape[0]} ')
# edges false need to be cut with the same shape as edges true, randomly
# check the section where we load edges_false_protein_drug_01 to see this
edges_all_false = edges_all_false.sample(n=edges_all.shape[0], random_state=123)
logging.info(f'# edges false now : {edges_all_false.shape[0]} ')
# this is how they "balance" the dataset, ramdomly chosing negatives

# save
edges_all.to_csv(os.path.join(path, f'{folder}/edges_all_(0, 1).txt'), sep= " ", header=None, index=None)
edges_all_false.to_csv(os.path.join(path, f'{folder}/edges_all_false(0, 1).txt'), sep= " ", header=None, index=None)


# now we need to create the folowing ("splitted" files)
# those are:
#   # 'train_index_(0, 1)i.txt'
#   # 'test_index_(0, 1)0.txt'
#   # 'index_test_false(0, 1)0.txt'

# we create the positive splits with

all_edge_idx = list(range(edges_all.shape[0]))
np.random.seed(123)
np.random.shuffle(all_edge_idx) # ? 

kf = KFold(n_splits=10, shuffle=True, random_state=123)
# Then we use the 10-flod cross-validation method to evaluate our mode
# we first generate the train set and test set.
i =0
for pos_train_idx, pos_test_idx in kf.split(all_edge_idx):
    # pos_train_idx is train set, pos_test_idx is test set
    ###### TRAIN  INDEX ( 0, 1 )
    train_edges = edges_all.iloc[pos_train_idx]  # These are  'train_index_(0, 1)0.txt'
    train_edges.to_csv(os.path.join(path, f'{folder}/train_index_(0, 1){i}.txt'), sep= " ", header=None, index=None)
    #
    ###### TEST  INDEX ( 0, 1 )
    test_edges = edges_all.iloc[pos_test_idx] 
    test_edges.to_csv(os.path.join(path, f'{folder}/test_index_(0, 1){i}.txt'), sep= " ", header=None, index=None)
    i+=1
# 


### negative test splits?
# we dont know exactly how they did
# so we repeat the process with the negative matrix a
# and only sabe those for test (only txt that the model takes)

all_edge_idx_false = list(range(edges_all_false.shape[0]))
np.random.seed(123)
np.random.shuffle(all_edge_idx_false) # ? 

# define kf again but actually not needed
# gonna be the same as the random state is 123
kf = KFold(n_splits=10, shuffle=True, random_state=123)

i =0
for pos_neg_train_idx, pos_neg_test_idx in kf.split(all_edge_idx_false):
    # They dont use the train neg edges 
    #train_neg_edges = edges_all_false.iloc[pos_neg_train_idx]  # The
    #print(train_neg_edges.shape)
    test_neg_edges = edges_all_false.iloc[pos_neg_test_idx]  # The
    print(test_neg_edges.shape)
    test_neg_edges.to_csv(os.path.join(path, f'{folder}/index_test_false(0, 1){i}.txt'), sep= " ", header=None, index=None)
    i+=1

# test
### checks 
# working with protein-drug
# 0 == PROTEIN
# 1 == DRUG
logging.info('-'*30)
edges_protein_drug_01 = pd.read_csv('luo/oneTooneIndex/edges_all_(0, 1).txt', sep=' ', header=None)
logging.info(f'edges_all is the same as the available file edges_all_(0,1): {(edges_all == edges_protein_drug_01).all().all()}')
# for false we need to check the file

edges_false_protein_drug_01 = pd.read_csv('luo/oneTooneIndex/edges_all_false(0, 1).txt', sep=' ', header=None)
edges_false_protein_drug_01.shape