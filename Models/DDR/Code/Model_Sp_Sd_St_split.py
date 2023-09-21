import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from itertools import product
from collections import defaultdict as dd
import random as r
import functools as f
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from random import randint
from functools import reduce as red

#Here we will define the model splits.
#We will have 4 type of splits, 3 of them have been already mentioned in the literature.
#We will be following DDR comparison scheme:
    # 5-repeats of 90%-10% 10-fold CV

### Sp, corresponds to the situation when there are
#DTI in the training data for such drugs or target proteins 
### Sd, corresponds to the situation when there are
#not DTI in the training data for some drugs
### St, corresponds to the situation when there are
#not DTI in the training data for some proteins
def generate_splits(DTIs, mode='Sp', subsampling=True, foldnum=10, cvopt=True, RMSD_dict = None):

    #For this split, we will use a regular Kfold split
    #Create Drug and Protein sets
    Drug_set = set(DTIs['Drug'].values)
    Drug_L = len(Drug_set)
    Drug_index = list(range(len(Drug_set)))

    Prot_set = set(DTIs['Protein'].values)
    Prot_L = len(Prot_set)
    Prot_index = list(range(Prot_L))

    #Create dictionaries for storing the positions
    Drug_dd = dict(zip(sorted(Drug_set),Drug_index))
    Drug_inv_dd = {v: k for k, v in Drug_dd.items()}
    Prot_dd = dict(zip(sorted(Prot_set),Prot_index))
    Prot_inv_dd = {v: k for k, v in Prot_dd.items()}

    #Create Drug-to-Proteins dict
    #DrugProteinDD, ProteinDrugDD = DrugToProteinDict(DTIs)
    
    def DrugToProteinDict(DTIs):

        #init drug to protein dd
        dpdd = dd(list)
        pddd = dd(list)

        #fill the dict
        for pair in DTIs.values.tolist():
            dpdd[pair[0]] += [pair[1]]

        for pair in DTIs.values.tolist():
            dpdd[pair[1]] += [pair[0]]

        return dpdd, pddd

    def get_interactions_dict(DTIs, seed, subsampling, swap = False, swap_dict = None, RMSD_dict = None):

        def get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict, Prot_inv_dd):

            def unify_genes(allItems, maxSample, negprots):

                
                #first sort
                sortAllItems= sorted(allItems, key = lambda x: x[1])

                #init Max
                uniqueSortedItems = []
                ref = sortAllItems[0][0]

                #remove duplicities
                for tupla in sortAllItems:
                    gen = tupla[0]

                    if gen != ref:
                        ref = gen
                        if gen in negprots:
                            uniqueSortedItems.append(gen)
                            if len(uniqueSortedItems) >= maxSample:
                                return uniqueSortedItems
                    
            #define maximum amount of genes to be sampled 
            maxSample = min(len(neg_element),len(pos_element))

            #get all proteins
            prots = [Prot_inv_dd[protid] for protid in pos_element]
            negprots =[Prot_inv_dd[protid] for protid in neg_element]

            #get all items
            allItems = []

            #concat
            for prot in prots:
                if prot in RMSD_dict:
                        allItems += RMSD_dict[prot].items()

            return unify_genes(allItems, maxSample, negprots)

        def get_drugs_for_targets_RMSD(pos_element, neg_element, RMSD_dict, Drug_inv_dd, DrugProteinDD, ProteinDrugDD):

            def unify_genes(allItems, maxSample, negdrugs, ProteinDrugDD):

                #first sort
                sortAllItems= sorted(allItems, key = lambda x: x[1])

                #init Max
                uniqueSortedItems = []
                ref = sortAllItems[0][0]

                #remove duplicities
                for tupla in sortAllItems:
                    gen = tupla[0]

                    if gen != ref:
                        ref = gen
                        if gen in negdrugs:
                            uniqueSortedItems.append(gen)
                            if len(uniqueSortedItems) >= maxSample:
                                return uniqueSortedItems

            #define maximum amount of genes to be sampled 
            maxSample = min(len(neg_element),len(pos_element))

            #get all proteins
            drugs = [Drug_inv_dd[protid] for protid in pos_element]
            negdrugs =[Drug_inv_dd[protid] for protid in neg_element]

            #get the associated proteins
            assocprots = set(red(lambda x,y : x+y, [DrugProteinDD[drug] for drug in drugs]))

            #get all items
            allItems = []

            #concat
            for prot in assocprots:
                if prot in RMSD_dict:
                        allItems += RMSD_dict[prot].items()

            return unify_genes(allItems,maxSample,negdrugs,ProteinDrugDD)

        #init default dict (list)
        interactions_dd = dd(list)

        #get prng
        prng = np.random.RandomState(seed)

        #get positives
        for d,p in zip(DTIs['Drug'],DTIs['Protein']):
            if swap:
                interactions_dd[Prot_dd[p]].append((Drug_dd[d],Prot_dd[p],1))
            else:
                interactions_dd[Drug_dd[d]].append((Drug_dd[d],Prot_dd[p],1))
               
        #add negatives (subsample to have 50%-50%)
        #go through all drugs/proteins
        for i, elementid in enumerate(interactions_dd):
            
            #print(f"elementid {elementid} in i {i}")
            #subsample from the negatives and add it to interactions dictionary
            #drugs if swap = False | proteins if swap = True
            if swap:
                pos_element = list(map(lambda x: x[0], interactions_dd[elementid])) #x[0] = drug
                neg_element = set(range(Drug_L)).difference((set(pos_element).union(swap_dict[elementid])))
            else:
                pos_element = list(map(lambda x: x[1], interactions_dd[elementid])) #x[1] = prot
                neg_element = set(range(Prot_L)).difference(set(pos_element))
            #print(f"Positive element {len(pos_element)}, negative elements {len(neg_element)}")

            if subsampling:
                #check if we can subsample all
                if RMSD_dict is None:
                    neg_sampled_element = r.sample(neg_element,min(len(neg_element),len(pos_element))) #50%-50% (modify if different proportions are desired)
                else:
                    if swap:
                        neg_sampled_element = get_drugs_for_targets_RMSD(pos_element, neg_element, RMSD_dict, Prot_inv_dd)
                    else:
                        neg_sampled_element = get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict, Prot_inv_dd)
            else:
                neg_sampled_element = neg_element #get all negatives

            #generate the negatives
            if swap:
                neg_sampled_element = [(drugid,elementid,0) for drugid in neg_sampled_element] #elementid is protid
            else:
                neg_sampled_element = [(elementid,protid,0) for protid in neg_sampled_element] #elementid is drugid
                
            #append
            interactions_dd[elementid] += neg_sampled_element
            #shuffle
            prng.shuffle(interactions_dd[elementid])

        return interactions_dd

    def set_to_matrix(set_interactions):
        return np.array(list(set_interactions),dtype=object)

    def names_from_edges(train_edges,test_edges):

        #generate names matrix from positive edges
        train_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p]) for d,p in train_edges]

        #same with negatives
        test_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p]) for d,p in test_edges]

        return train_names_matrix, test_names_matrix

    def optimize_folds(cv_distribution, pos_neg_interactions, pos_neg_interactions_dd):
        #here we will distribute drugs according to Sd
        #maximizing the number of different drugs we have in each fold

        #first we compute length of each drug
        drugs_L_tuple = sorted([(drug,len(pos_neg_interactions_dd[drug])) for drug in pos_neg_interactions_dd], key= lambda x: x[1])

        #compute elements per fold
        elements_per_fold = len(pos_neg_interactions)//foldnum

        i = 0
        while True:
            drugs_tuple = drugs_L_tuple.pop(0)
            #length
            fold_length = drugs_tuple[1]
            #check if its empty enough to include that drug
            if len(cv_distribution[i%foldnum]) + fold_length <= elements_per_fold:
                cv_distribution[i%foldnum]+=pos_neg_interactions_dd[drugs_tuple[0]]
            else:
                for element in pos_neg_interactions_dd[drugs_tuple[0]]:
                    cv_distribution[i%foldnum].append(element)
            if not len(drugs_L_tuple):
                break
            i+=1

        return cv_distribution

    #init seed cv list
    seed_cv_list = []

    for seed in tqdm([7183,556,2,81,145], desc='Performing 10-CV fold for each seed'):

        # check if subsampling
        if subsampling:

            ##drugs
            drug_interactions_dd = get_interactions_dict(DTIs, seed, subsampling=True, swap=False)
            # append all interactions
            pos_neg_interactions_drugs = f.reduce(lambda a,b: a+b, drug_interactions_dd.values())

            #generate the inverse drug_interaction_dd to avoid duplicates when 
            #merging drug/proteins interactions_dd
            def inverse_interactions_dd(drug_interactions_dd):
                interactions_dd = dd(set)

                for drugid in drug_interactions_dd:
                    for element in drug_interactions_dd[drugid]:
                        #add entry (protein) = {drugs that have already been used}
                        interactions_dd[element[1]].add(drugid)
            
                return interactions_dd

            #compute inverse dict
            inv_drug_interactions_dd = inverse_interactions_dd(drug_interactions_dd)

            ##proteins
            prot_interactions_dd = get_interactions_dict(DTIs, seed, subsampling=True, swap=True, swap_dict = inv_drug_interactions_dd)
            # append all interactions
            pos_neg_interactions_proteins = f.reduce(lambda a,b: a+b, prot_interactions_dd.values())

            def merge_interactions(drugs_interactions,prots_interactions):

                #merge
                interactions = drugs_interactions + prots_interactions
                #remove possible duplicities
                unique_interactions = []
                #(use a dict for checking)
                interactions_d = {}
                for element in interactions:
                    if element not in interactions_d:
                        interactions_d[element] = None
                        unique_interactions.append(element)

                return unique_interactions

            #merge
            pos_neg_interactions = merge_interactions(pos_neg_interactions_drugs,pos_neg_interactions_proteins)
           
        else:
            drug_interactions_dd = get_interactions_dict(DTIs, seed, subsampling=False)
            # append all interactions
            pos_neg_interactions = f.reduce(lambda a,b: a+b, drug_interactions_dd.values())
        
        #check % of positives/negatives
        pos_percentage = sum(list(map(lambda x: x[2],pos_neg_interactions)))/len(pos_neg_interactions)
        print(f"Positives -> {round(pos_percentage,2)*100}%, Negatives -> {round(1-pos_percentage,2)*100} %")

        #init list to distribute edges in a Sp way.
        cv_distribution = [[] for _ in range(foldnum)]

        if mode == 'Sp':

            #both targets and drugs have been already seen during the training, but this exact DTI is new.
            for i,interaction in enumerate(tqdm(pos_neg_interactions,desc='Distributing interactions')):
                #get positives for that drug
                cv_distribution[i%foldnum].append(interaction)

        elif mode == 'Sd':

            #drugs are new, targets have been seen during the training
            pos_neg_interactions_dd = dd(list)
            for interaction in pos_neg_interactions:
                pos_neg_interactions_dd[interaction[0]].append(interaction)

            cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions, pos_neg_interactions_dd)

        elif mode == 'St':

            pos_neg_interactions_dd = dd(list)
            #prots are new, drugs have been seen during the training
            for interaction in pos_neg_interactions:
                pos_neg_interactions_dd[interaction[1]].append(interaction)

            cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions, pos_neg_interactions_dd)


        #generate the interaction matrix
        pos_neg_matrix = set_to_matrix(f.reduce(lambda a,b: a+b, cv_distribution))


        ## MIKEL OPTION
        if not cvopt:

            X_train, X_val_test = train_test_split(pos_neg_matrix, test_size=0.3, random_state=None)
            X_val, X_test = train_test_split(X_val_test, test_size=0.66, random_state=None)

            def get_arkar():
                train_edges_pos, test_edges_pos = X_train[X_train[:,2] == 1,:-1], X_train[X_train[:,2] == 1,:-1]   
                names_train_pos, names_test_pos = names_from_edges(train_edges_pos,test_edges_pos)

            #1 seed -> [[pos_train,neg_train],[pos_val,neg_val],[pos_test,neg_test]]
            #2 seed -> ...
            #..
            #5 seed -> ...

        else:

            #define both KFold to maintain proportions (positive and negatives)
            KFoldobj = KFold(n_splits = foldnum)

            #init the cv list
            cv_list = []

            for train_index, test_index in KFoldobj.split(pos_neg_matrix):
            
                #get train,test
                train_edges, test_edges = pos_neg_matrix[train_index],pos_neg_matrix[test_index]
                
                #--positives--
                train_edges_pos, test_edges_pos = train_edges[train_edges[:,2] == 1,:-1], test_edges[test_edges[:,2] == 1,:-1]
                ##create names matrix from edges list
                names_train_pos, names_test_pos = names_from_edges(train_edges_pos,test_edges_pos)

                #--negatives--
                train_edges_neg, test_edges_neg = train_edges[train_edges[:,2] == 0,:-1], test_edges[test_edges[:,2] == 0,:-1]
                
                ##create names matrix from edges list
                names_train_neg, names_test_neg = names_from_edges(train_edges_neg,test_edges_neg)

                #print(f"Train pos {len(names_train_pos)}, Train neg {len(names_train_neg)}, Test pos {len(names_test_pos)}, Test neg {len(names_test_neg)}")

                #add each fold
                cv_list.append((names_train_pos,names_train_neg,names_test_pos,names_test_neg))

            #add each group of folds for each seed
            seed_cv_list.append(cv_list)

    return seed_cv_list

def genRMSDdict():

    fpath = '/mnt/md0/data/jfuente/DTI/Input4Models/Docking/Results/RMSD_full_matrix.pkl'

    #load RMSD object
    RMSD = pd.read_pickle(fpath)

    #init dict and names
    RMSDdict = dict()
    names = RMSD.index.tolist()

    # convert to numpy
    RMSD = RMSD.to_numpy()

    #fill it

    for i,gen in enumerate(tqdm(names)):

        #get genes and names
        rmsd_i = RMSD[i,0:i].tolist() + RMSD[i,i+1:].tolist()
        names_i = names[0:i] + names[i+1:]

        #add the entry
        RMSDdict[gen] = dict(zip(names_i,rmsd_i))


    return RMSDdict

#Lets use Yamanishi NR as an example
##Load dataset
fpath = os.path.join(os.getcwd(),'DB','Data','Yamanashi_et_al_GoldStandard','IC','interactions','ic_admat_dgc_mat_2_line.txt')
#fpath = os.path.join(os.getcwd(),'DB','Data','Davis_et_al','tdc_package_preprocessing','DAVIS_et_al_2line.tsv')
#fpath = os.path.join(os.getcwd(),'DB','Data','BIOSNAP','ChG-Miner_miner-chem-gene','ChG-Miner_miner-chem-gene.tsv')
DTIs = pd.read_csv(fpath,sep='\t')
DTIs.columns = ['Drug','Protein']

#check splits
def check_splits(splits, verbose=False, foldnum=10):

    def check_condition():
        drugs_counter, prots_counter = 0,0
        for seed in range(len(splits)):
            #print(f"Seed {seed}")
            for fold in range(len(splits[seed])):
                #TRAIN
                drugs_train = set(map(lambda x: x[0], splits[seed][fold][0])).union(map(lambda x: x[0], splits[seed][fold][1]))
                prots_train = set(map(lambda x: x[1], splits[seed][fold][0])).union(map(lambda x: x[1], splits[seed][fold][1]))
                #TEST
                drugs_test = set(map(lambda x: x[0], splits[seed][fold][2])).union(map(lambda x: x[0], splits[seed][fold][3]))
                prots_test = set(map(lambda x: x[1], splits[seed][fold][2])).union(map(lambda x: x[1], splits[seed][fold][3]))

                if drugs_test.difference(drugs_train):
                    drugs_counter +=1
                    if verbose:
                        print(f"Seed {seed}, Fold {fold} does not accomplish drugs")

                if prots_test.difference(prots_train):
                    prots_counter+=1
                    if verbose:
                        print(f"Seed {seed}, Fold {fold} does not accomplish proteins")


        if drugs_counter > 0:
            if drugs_counter == (len(splits)*foldnum):
                print("Sd split configuration accomplished!")

            else:
                print("Sd split configuration not exactly accomplished!")

        if prots_counter > 0:
            
            if prots_counter == (len(splits)*foldnum):
                print("St split configuration accomplished!")
            else:
                print("St split configuration not exactly accomplished!")

        if not prots_counter and not drugs_counter:
            print('Sp split configuration accomplished!')

    def print_proportions(verbose=True):
        for seed in range(len(splits)):
            if verbose:
                print(f"Len is {len(splits[seed])}")
            for fold in range(len(splits[seed])):
                if verbose:
                    print(f"Train positives len {len(splits[seed][fold][0])}")
                    print(f"Train negatives len {len(splits[seed][fold][1])}")
                    print(f"Test positives len {len(splits[seed][fold][2])}")
                    print(f"Test negatives len {len(splits[seed][fold][3])}")
                    print("--------------- END OF FOLD ---------------")
            if verbose:
                print("----------------- END OF SEED -----------------")

    print('Checking conditions per fold')
    #check condition
    check_condition()

    # if verbose:
    #     print('Printing proportions')
    # #print proportions
    # print_proportions(verbose=verbose)

    return 

#check distribution
def print_cv_distribution(DTIs, cv_distribution):

    print(f'Number of drugs originally -> {len(set(DTIs.Drug.values))}')
    print(f'Number of prots originally -> {len(set(DTIs.Protein.values))}')

    element_distribution = f.reduce(lambda a,b: a+b, cv_distribution)
    drugs = set(list(map(lambda x: x[0], element_distribution)))
    prots = set(list(map(lambda x: x[1], element_distribution)))

    print(f'Number of drugs in CV -> {len(drugs)}')
    print(f'Number of prots in CV-> {len(prots)}')

    for drug in drugs:
        proteins_with_drugs = [element[1] for element in element_distribution if element[0] == drug]
        print(f"Number of proteins per drug {drug} -> {len(proteins_with_drugs)}")

    for prot in prots:
        drugs_with_proteins = [element[0] for element in element_distribution if element[1] == prot]
        print(f"Number of drugs per protein {prot} -> {len(drugs_with_proteins)}")

    ##~
    return 

#generate RMSD dictionary
RMSD_dict = genRMSDdict()

# ------------------------------------------------------- Sp ------------------------------------------------------------------ #
##Get 5-seed 10-fold CV Sp (all nodes are seeing during the training)
sp_splits = generate_splits(DTIs, mode= 'Sp', subsampling=True, foldnum=10)

##check sp split
check_splits(sp_splits,verbose=False) 

##check distribution
print_cv_distribution(DTIs, sp_splits[0])

# ------------------------------------------------------- Sd ------------------------------------------------------------------- #
##Get 5-seed 10-fold CV Sd (some drugs are not seeing during the training)
sd_splits = generate_splits(DTIs, mode= 'Sd', subsampling=True, foldnum=10)

##check splits
check_splits(sd_splits,verbose=False) 

##check distributions
#print_cv_distribution(DTIs,sd_splits)
# ------------------------------------------------------- St ------------------------------------------------------------------- #
##Get 5-seed 10-fold CV Sd (some targets are not seeing during the training)
st_splits = generate_splits(DTIs, mode= 'St', subsampling=True, foldnum=10)

##check splits
check_splits(st_splits,verbose=False) 