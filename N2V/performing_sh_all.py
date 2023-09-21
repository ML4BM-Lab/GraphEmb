import pandas as pd
import os
import tensorflow as tf
import sys
import importlib
# import models
from tensorflow.keras import layers
import tensorflow_addons as tfa
import numpy as np
import random
import itertools
import sklearn
from tqdm import tqdm
import pickle
from sklearn.model_selection import train_test_split

#sys.path.append(os.getcwd()+'/models/gold_standard_classificator/')
#importlib.reload(models)

db_dict = './DBs_dict/'
emb_p = './DBs_embeddings/'

"""
Yamanishi-E
BIOSNAP
"""

#read dataframe
n2v_results = pd.read_csv(os.path.join('panel','figure1','n2v_nn_results.tsv'), sep='\t', index_col = 0)


all_models = ['Yamanishi/E','Yamanishi/NR','Yamanishi/GPCR','Yamanishi/IC',
              'BIOSNAP','BindingDB','DrugBank','Davis_et_al']

def train_model(X_train, y_train, X_validation, y_validation, comb):

    emb_dim = X_train.shape[1]

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
    batch_size = round(X_train.shape[0]*batch_size)

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
    
    #compute and append to list
    train_auprc = sklearn.metrics.average_precision_score(y_train, y_score_train)
    val_auprc = sklearn.metrics.average_precision_score(y_validation, y_score_validation)
    ## AUC
    #Compute AUC and validation AUC
    train_auc, val_auc = history.history['auc'][-1], history.history['val_auc'][-1]

    return train_auc, val_auc, train_auprc, val_auprc

def retrieve_data(db_chosen, best_run_title):

    if 'Yamanishi' in best_run_title:
        i = 1
    elif 'Davis' in best_run_title:
        i = 2
    else:
        i = 0

    dims = best_run_title.split('_')[1+i]
    pvar = best_run_title.split('_')[2+i]
    lvar = best_run_title.split('_')[3+i]

    db_chosen_edges = db_chosen.replace('/','_') if 'Yamanishi' in db_chosen else db_chosen

    #load edges
    edges = pd.read_csv(f'{db_dict}{db_chosen_edges}_edges.tsv',sep='\t',header=None)
    drugs = edges.iloc[:,0].values
    genes = edges.iloc[:,1].values

    #generate the positive connections
    positive_connections = set([tuple(x) for x in zip(drugs, genes)])

    #generate the negative connections (we always have GENE-DRUG)
    negative_connections = random.sample(set(itertools.product(drugs, genes)).difference(positive_connections), len(positive_connections))

    #build total connections
    total_connections = list(positive_connections)+list(negative_connections)

    #build the y
    y = np.array([1]*len(positive_connections)+[0]*len(negative_connections))


    emb_f_p = f'{db_chosen_edges}_{dims}_{pvar}_{lvar}_embedding.emb' 
    emb_f = pd.read_csv(f'{emb_p}{db_chosen}/{emb_f_p}',skiprows = [0],sep=' ',header=None)

    #build a dict containing (ids,embeddings)
    emb_f_dict = dict(zip(emb_f.iloc[:,0],emb_f.iloc[:,1:].values.tolist()))

    #build the (X)
    X = np.array([emb_f_dict[edge[0]]+emb_f_dict[edge[1]] for edge in total_connections])

    return X,y

all_results = []

if __name__ == "__main__":
    
    for model_train in tqdm(all_models, desc='looping through all the models'):
    
        model_parsed = model_train.replace('/','_') if 'Yamanishi' in model_train else model_train
        model_parsed = 'Davis' if 'Davis' in model_parsed else model_parsed

        #get best result for model according to validation AUC
        model_best = n2v_results[n2v_results['Network'].values == model_parsed]['Val AUC'].sort_values(ascending=False).index[0]

        #read info
        batch_size_v = float(model_best.split('_')[-5])
        epochs_v = int(model_best.split('_')[-4])
        model_type = '_'.join(model_best.split('_')[-3:-1])
        loss_type = model_best.split('_')[-1].split('.')[0]

        #define the train group of parameters
        comb = (batch_size_v, epochs_v, model_type, loss_type)

        #
        X_train, y_train = retrieve_data(model_train, model_best)

        for model_val in all_models:
            
            if model_val == model_train:
                X_train, X_validation, y_train, y_validation = train_test_split(X_train, y_train, test_size = 0.3)
            else:
                X_validation, y_validation = retrieve_data(model_val, model_best)

            #load the model
            train_auc, val_auc, train_auprc, val_auprc = train_model(X_train, y_train, X_validation, y_validation, comb)

            all_results.append([model_train,model_val,(train_auc, val_auc, train_auprc, val_auprc)])

            print(f"""When training on {model_train} and testing on {model_val}, results
                    were: train AUC {train_auc}, validation AUC {val_auc},
                    train AUPRC {train_auprc}, validation AUPRC {val_auprc}""")

with open(os.path.join('figure2','all_results.pickle'), 'wb') as handle:
    pickle.dump(all_results, handle, protocol=pickle.HIGHEST_PROTOCOL)