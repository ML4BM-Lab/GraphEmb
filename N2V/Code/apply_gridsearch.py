import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import json
import pandas as pd
import itertools
import numpy as np
from imblearn.over_sampling import SMOTE, RandomOverSampler
from collections import Counter
import importlib
import sys
import pickle
import sklearn.metrics
from tqdm import tqdm
import random
import tensorflow as tf
from tensorflow.keras import layers
import tensorflow_addons as tfa
import os
import json
import pandas as pd
import itertools
import numpy as np
from collections import Counter
import copy
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import binarize
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

db_dict = 'Data/DBs_dict/'
emb_p = 'Data/DBs_embeddings/'

#Define the parameters we are going to be varying for the grid
batch_size_v = [1/16, 1/64]
epochs_v = [2, 5, 10, 20, 50]
model_type = ['Type_1', 'Type_2', 'Type_3']
loss_type = ['BCE', 'Focal']

def train_model(X,y,comb):

    emb_dim = X.shape[1]

    batch_size, epochs,model_type,loss_type = comb

    #Here we are defining some models (4), all based in the emb_dim
    #Model
    if model_type=='Type_1':
        drug_gene_model = tf.keras.Sequential([
            layers.Dense(units=round(emb_dim/2),activation='relu',input_shape=(emb_dim,)),
            layers.Dense(units=1,activation='sigmoid')
        ])
    elif model_type=='Type_2':
        drug_gene_model = tf.keras.Sequential([
            layers.Dense(units=round(emb_dim),activation='relu',input_shape=(emb_dim,)),
            layers.Dense(units=1,activation='sigmoid')
        ])
    elif model_type=='Type_3':
        drug_gene_model = tf.keras.Sequential([
            layers.Dense(units=round(emb_dim/2),activation='relu',input_shape=(emb_dim,)),
            layers.Dense(units=round(emb_dim/4),activation='relu'),
            layers.Dense(units=1,activation='sigmoid')
        ])
    elif model_type=='Type_4':
        drug_gene_model = tf.keras.Sequential([
            layers.Dense(units=emb_dim,activation='relu',input_shape=(emb_dim,)),
            layers.Dense(units=round(emb_dim/2),activation='relu'),
            layers.Dense(units=1,activation='sigmoid')
        ])
    
    #LOSS TYPE
    #If we use L2 regularization, we use BCE instead of FocalCE
    if loss_type=='BCE':
        drug_gene_model.compile(optimizer='adam',
        loss = tf.keras.losses.BinaryCrossentropy(from_logits=False),
        metrics=['binary_crossentropy','accuracy','TrueNegatives',
                'TruePositives','FalsePositives','FalseNegatives',tfa.metrics.F1Score(num_classes=1),'AUC'])
    elif loss_type=='Focal': #Use FocalCE otherwise
        drug_gene_model.compile(optimizer='adam',
        loss= tfa.losses.SigmoidFocalCrossEntropy(from_logits=False),
        metrics=['binary_crossentropy','accuracy','TrueNegatives',
                'TruePositives','FalsePositives','FalseNegatives',tfa.metrics.F1Score(num_classes=1),'AUC'])

    #Process batch_size
    batch_size = round(X.shape[0]*batch_size)

    ## TRAIN || VALIDATION || TEST SPLITS
    train_ratio = 0.75
    validation_ratio = 0.15
    test_ratio = 0.1

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 1 - train_ratio, random_state=42)
    X_validation, X_test, y_validation, y_test = train_test_split(X_test, y_test, test_size = test_ratio/(test_ratio + validation_ratio), random_state=1)

    #TRAIN
    history = drug_gene_model.fit(
                    x=X_train,
                    y=y_train,
                    epochs=epochs,
                    batch_size=batch_size,
                    validation_data = (X_validation,y_validation),
                    shuffle=True,
                    verbose=0
                    )
    

    ## AUPRC
    #first compute y_pred
    y_score_train = drug_gene_model.predict(X_train)
    y_score_validation = drug_gene_model.predict(X_validation)
    y_score_test = drug_gene_model.predict(X_test)
    
    #compute and append to list
    train_auprc = sklearn.metrics.average_precision_score(y_train, y_score_train)
    val_auprc = sklearn.metrics.average_precision_score(y_validation, y_score_validation)
    test_auprc = sklearn.metrics.average_precision_score(y_test, y_score_test)

    ## AUC
    #Compute AUC and validation AUC
    train_auc, val_auc = history.history['auc'][-1], history.history['val_auc'][-1]
    test_auc = roc_auc_score(y_test, y_score_test)


    return train_auc, val_auc, test_auc, train_auprc, val_auprc, test_auprc

#choose the database we want to work with
db_chosen = sys.argv[1]

if __name__ == "__main__":

    db_chosen_edges = db_chosen.replace('/','_') if 'Yamanishi' in db_chosen else db_chosen

    #load edges
    edges = pd.read_csv(f'{db_dict}{db_chosen_edges}_edges.tsv',sep='\t',header=None)
    drugs = edges.iloc[:,0].values
    genes = edges.iloc[:,1].values

    """
    NOT NECESSARY FOR NOW
    #build drug dicts
    drug_dict = pd.read_csv(f'{db_dict}{db_chosen}drug_dict.tsv',sep='\t',header=None)

    #build prot dicts
    prot_dict = pd.read_csv(f'{db_dict}{db_chosen}prot_dict.tsv',sep='\t',header=None)
    """

    #generate the positive connections
    positive_connections = set([tuple(x) for x in zip(drugs, genes)])

    #generate the negative connections (we always have GENE-DRUG)
    negative_connections = random.sample(set(itertools.product(drugs, genes)).difference(positive_connections), len(positive_connections))

    #build total connections
    total_connections = list(positive_connections)+list(negative_connections)

    #build the y
    y = np.array([1]*len(positive_connections)+[0]*len(negative_connections))

    for emb_f_p in tqdm(os.listdir(emb_p+db_chosen+'/'), desc = f'Going through embeddings from network {db_chosen}'):
        
        emb_f = pd.read_csv(f'{emb_p}{db_chosen}/{emb_f_p}',skiprows = [0],sep=' ',header=None)

        #build a dict containing (ids,embeddings)
        emb_f_dict = dict(zip(emb_f.iloc[:,0],emb_f.iloc[:,1:].values.tolist()))

        #build the (X)
        X = np.array([emb_f_dict[edge[0]]+emb_f_dict[edge[1]] for edge in total_connections])

        for comb in tqdm(list(itertools.product(batch_size_v,epochs_v,model_type,loss_type)), desc = 'Hyperparameter Tunning Grid'):

            #get the identifier (old name + new tags)
            ID = emb_f_p[:-4]+'_'+'_'.join(map(str,comb))

            file_path = os.getcwd()+f'/Results/{ID}.pickle'

            #check if file exists
            if os.path.isfile(file_path):
                continue

            #train the model 
            train_auc, val_auc, test_auc, train_auprc, val_auprc, test_auprc = train_model(X,y,comb)
            summary_v = [train_auc, val_auc, test_auc, train_auprc, val_auprc, test_auprc]

            #save dict
            with open(file_path, 'wb') as handle:
                pickle.dump(summary_v, handle, protocol=pickle.HIGHEST_PROTOCOL)
