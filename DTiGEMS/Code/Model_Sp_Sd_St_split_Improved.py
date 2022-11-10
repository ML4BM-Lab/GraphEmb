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
from collections import Counter 
import itertools
import math as m

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
def generate_splits(DTIs, mode='Sp', subsampling=True, foldnum=10, negative_to_positive_ratio = 1,  
                    cvopt=True, cvopt_no_names = False, ttv_names = True, 
                    train_val_test_percentage = (0.7, 0.1, 0.2), 
                    RMSD_dict_opt = False, RMSD_threshold = 6, only_distribution = False, 
                    include_diagonal_RMSD = False, n_seeds = 5):

    def genRMSDdict(genes):

        def filterRMSD(RMSD, DTIs, names):

            #get DTIs proteins
            DTIsprots = DTIs['Protein'].values

            #get the index
            ind = [i for i,x in enumerate(names) if x in DTIsprots]

            return RMSD[ind, :][:,ind], list(np.array(names)[ind])

        fpath = '/mnt/md0/data/jfuente/DTI/Input4Models/Docking/Results/RMSD_full_matrix.pkl'

        #load RMSD object
        RMSD = pd.read_pickle(fpath)

        #init dict and names
        RMSDdict = dict()
        names = RMSD.index.tolist()

        # convert to numpy
        RMSD = RMSD.to_numpy()

        #Make sure every entry in RMSD matrix is in DTIs
        RMSD, names = filterRMSD(RMSD, DTIs, names)

        #fill it
        for i,gen in enumerate(tqdm(names)):

            if gen not in genes:
                continue

            if include_diagonal_RMSD:
                #get genes and names
                rmsd_i = RMSD[i,:].tolist()
                names_i = names
            else:
                #get genes and names
                rmsd_i = RMSD[i,0:i].tolist() + RMSD[i,i+1:].tolist()
                names_i = names[0:i] + names[i+1:]

            #add the entry
            RMSDdict[gen] = dict(zip(names_i,rmsd_i))

        return RMSDdict

    def FilterDTIs(DTIs,RMSD_dict):

        print(f"Original shape: {len(set(DTIs['Drug']))} drugs x {len(set(DTIs['Protein']))} proteins")

        fdrug = []
        ftarget = []
        for drug, target in zip(DTIs['Drug'], DTIs['Protein']):
            if target in RMSD_dict:
                fdrug.append(drug)
                ftarget.append(target)

        fDTIs = pd.DataFrame([fdrug,ftarget]).T
        fDTIs.columns=['Drug','Protein']

        print(f"Filtered shape: {len(set(fDTIs['Drug']))} drugs x {len(set(fDTIs['Protein']))} proteins")

        return fDTIs

    if RMSD_dict_opt:

        def get_positives_for_final_fold(DTIs):

            #Reset index just in case
            DTIs.reset_index(drop=True,inplace=True)

            #Build dict drug-to-proteins
            drug_to_prots = dd(list)
            tupla_index = {}
            final_pos_edges = []
            drop_prots = []
            
            for i,tupla in enumerate(zip(DTIs['Drug'], DTIs['Protein'])):
                drug_to_prots[tupla[0]].append(tupla[1])
                tupla_index[tupla] = i

            #Build also a protein counter, so when we move edges tp the final fold we are not dropping out the protein within that edge
            protein_d = dict(Counter(DTIs['Protein']))

            for drug in drug_to_prots:

                if len(drug_to_prots[drug]) > 2: 
                    
                    #compute protein multiplicity for the evaluated drug
                    protein_multiplicity = list(map(lambda x: protein_d[x], drug_to_prots[drug]))
                    if np.max(protein_multiplicity) > 5:
                        
                        #choose the protein with the most multiplicity
                        chosen_protein = drug_to_prots[drug][np.argmax(list(map(lambda x: protein_d[x], drug_to_prots[drug])))]

                        #update the multiplicity of that protein
                        protein_d[chosen_protein] -= 1

                        #append the chosen edge to drop it
                        final_pos_edges.append((drug,chosen_protein,1))
                        drop_prots.append(tupla_index[final_pos_edges[-1][:-1]])

            #Define the kept edges 
            kept_edges = list(set(range(DTIs.shape[0])).difference(drop_prots))

            DTIs_kept = DTIs.loc[kept_edges,:]
            DTIs_kept.reset_index(drop=True,inplace=True)
            DTIs_kept_l = list(zip(DTIs_kept['Drug'], DTIs_kept['Protein']))

            # Check there is no DTI 
            assert all([out_of_sample_edge[1] not in DTIs_kept_l for out_of_sample_edge in final_pos_edges])

            print(f"Shape after dropping out some positive proteins {len(set(DTIs_kept['Drug'].values))} x {len(set(DTIs_kept['Protein'].values))}")

            return final_pos_edges, DTIs_kept

        print("Applying RMSD sim matrix to perform subsampling!")
        init_genes = set(DTIs['Protein'].values)
        RMSD_dict = genRMSDdict(init_genes)
        DTIs = FilterDTIs(DTIs, RMSD_dict)

        ## Define a final fold to test the RMSD
        positive_final_fold, DTIs = get_positives_for_final_fold(DTIs)
        negative_final_fold = []

        
    if not cvopt:
        train_ratio = train_val_test_percentage[0]
        validation_ratio = train_val_test_percentage[1]
        test_ratio = train_val_test_percentage[2]
        foldnum = m.ceil(1/test_ratio)

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

    def get_interactions_dict(DTIs, seed, subsampling):

        def get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict):

            def unify_genes(allItems, negprots, sampnegprots):

                #first sort
                sortAllItems = sorted(allItems, key = lambda x: x[1])

                """
                Before applying the boxes' method
                ---------------------------------
                

                sortAllItems_thresholded = [prot_rmsd_tuple for prot_rmsd_tuple in sortAllItems if prot_rmsd_tuple[1] > RMSD_threshold]

                #Get the kept out items and select the last one (closest to the threshold)
                keptOutItems = sorted(set(sortAllItems).difference(sortAllItems_thresholded), key = lambda x: x[1])
                negative_final_fold.append((Drug_inv_dd[elementid], keptOutItems[-1][0],0))

                """
                """
                After applying the boxe's method
                """

                train_val_lth = 5
                train_val_uth = RMSD_threshold

                sortAllItems_thresholded = [prot_rmsd_tuple for prot_rmsd_tuple in sortAllItems if prot_rmsd_tuple[1] > train_val_lth and prot_rmsd_tuple[1] < train_val_uth]
                
                if not len(sortAllItems_thresholded):
                    return

                r.shuffle(sortAllItems_thresholded)

                ## final fold
                final_fold_lth = 2.5
                final_fold_uth = 5
                #Get the kept out items and select the last one (closest to the threshold)
                keptOutItems = [prot_rmsd_tuple for prot_rmsd_tuple in sortAllItems if prot_rmsd_tuple[1] > final_fold_lth and prot_rmsd_tuple[1] < final_fold_uth]
                if len(keptOutItems):
                    negative_final_fold.append((Drug_inv_dd[elementid], r.sample(keptOutItems,1)[0][0],0))
                

                #remove duplicities
                for tupla in sortAllItems_thresholded:
                    #print(tupla)
                    gen = tupla[0]
                    #print(f"RMSD {tupla[1]}")
                    #check the gen is in negative proteins and not in already chosen negative proteins
                    if Prot_dd[gen] not in sampnegprots and gen in negprots:
                        return gen
                    
            #define maximum amount of genes to be sampled 
            maxSample = min(len(neg_element), len(pos_element))

            #get all proteins
            prots = [Prot_inv_dd[protid] for protid in pos_element]
            negprots = [Prot_inv_dd[protid] for protid in neg_element]

            #get all items
            sampnegprots = []

            #concat
            for prot in prots:

                if maxSample == 0:
                    break

                if include_diagonal_RMSD:
                    sampled_negative = prot
                else:
                    sampled_negative = unify_genes(RMSD_dict[prot].items(), negprots, sampnegprots)

                if sampled_negative:
                    sampnegprots.append(Prot_dd[sampled_negative])

                maxSample -= 1

            return sampnegprots

        #init default dict (list)
        interactions_dd = dd(list)

        #get prng
        prng = np.random.RandomState(seed)

        #get positives
        for d,p in zip(DTIs['Drug'],DTIs['Protein']):
            interactions_dd[Drug_dd[d]].append((Drug_dd[d],Prot_dd[p],1))
               
        #add negatives (subsample to have 50%-50%)
        #go through all drugs/proteins
        for i, elementid in enumerate(tqdm(interactions_dd)):
            #print(f"elementid {elementid} in i {i}")
            #subsample from the negatives and add it to interactions dictionary
            #drugs if swap = False | proteins if swap = True
            pos_element = list(map(lambda x: x[1], interactions_dd[elementid])) #x[1] = prot
            neg_element = set(range(Prot_L)).difference(set(pos_element))
            #print(f"Positive element {len(pos_element)}, negative elements {len(neg_element)}")

            if subsampling:
                #check if we can subsample all
                if not RMSD_dict_opt:
                    neg_sampled_element = r.sample(neg_element, min(len(neg_element), negative_to_positive_ratio * len(pos_element))) #50%-50% (modify if different proportions are desired)
                else:
                    neg_sampled_element = get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict)
            else:
                neg_sampled_element = neg_element #get all negatives

            #generate the negatives
            neg_sampled_element = [(elementid,protid,0) for protid in neg_sampled_element] #elementid is drugid
                
            #append
            interactions_dd[elementid] += neg_sampled_element

            #shuffle
            prng.shuffle(interactions_dd[elementid])

        return interactions_dd

    def set_to_matrix(set_interactions):
        return np.array(list(set_interactions),dtype=object)

    def names_from_edges(train_edges,test_edges):

        if cvopt:
            #generate names matrix from positive edges
            train_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p]) for d,p in train_edges]
            #same with negatives
            test_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p]) for d,p in test_edges]
        else:
            #generate names matrix from positive edges
            train_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p], label) for d,p,label in train_edges]
            #same with negatives
            test_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p], label) for d,p,label in test_edges]

        return train_names_matrix, test_names_matrix

    def optimize_folds(cv_distribution, pos_neg_interactions_dd):
        #here we will distribute drugs according to Sd
        #maximizing the number of different drugs we have in each fold

        #first we compute length of each drug
        drugs_L_tuple = sorted([(drug,len(pos_neg_interactions_dd[drug])) for drug in pos_neg_interactions_dd], key= lambda x: x[1], reverse=True)

        i = 0
        while True:
            drugs_tuple = drugs_L_tuple.pop(0)

            #elements to add
            elems = pos_neg_interactions_dd[drugs_tuple[0]]

            #add to cv_distr
            cv_distribution[i%foldnum] += elems
            
            if not len(drugs_L_tuple):
                break

            i+=1

        return cv_distribution

    def Sp_reorder(cv_distribution):

        #We are going to divide DTIs in 2 groups, crucial and non-crucial
        #For a subsampling dataset to accomplish Sp requirement:

        #Drugs with +=2 targets > Drugs with == 1 target
        #(OR)
        #Proteins with +=2 drugs > Proteins with == 1 drugs

        ## Non-Crucial DTI (Di - Tj) requirement
        # Di is at least in 3 folds (2 folds would accomplish Sp requirement ✓, so if its in 3 then if moved will still accomplish)
        # Tj also accomplish this requirement

        ## Crucial DTI (Di - Tj) requirement
        # Tj is only present in 1 fold, so it needs to be swap with a Non-Crucial DTI

        def compute_cruciality(cv_distribution, pos_neg_interactions):

            def build_fold_dict(cv_distribution):

                def bfd_helper(mode='drug'):

                    if mode=='drug':
                        ind = 0
                    elif mode == 'prot':
                        ind = 1

                    #init default dict  (default int is 0)
                    fold_dd = dd(lambda: dd(int))

                    for i,fold in enumerate(cv_distribution):
                        fold_elements = set(list(map(lambda y: y[ind], fold)))
                        for element in fold_elements:
                            fold_dd[i][element] += 1
                        
                    return fold_dd

                protfold_dd = bfd_helper(mode='prot')
                drugfold_dd = bfd_helper(mode='drug')

                return drugfold_dd, protfold_dd

            def cc_helper(cv_distribution, fold_dd, elements):

                #init crucial and non-crucial lists
                crucial = []
                non_crucial = []

                for elm in elements:

                    prot_trues = sum([1 if fold_dd[i][elm] >= 1 else 0 for i,_ in enumerate(cv_distribution)])

                    if prot_trues >= 3:
                        non_crucial.append(elm)
                    elif prot_trues == 1:
                        crucial.append(elm)

                    # debugging
                    # if prot in [336, 682, 775, 1502]:
                    #     print(prot_trues)

                return crucial, non_crucial

            #get the proteins, as the drugs are well distributed by the way we have generated cv_distribution
            prots = set(list(map(lambda x: x[1], pos_neg_interactions)))
            drugs = set(list(map(lambda x: x[0], pos_neg_interactions)))

            #retrieve prot fold dd
            drugfold_dd, protfold_dd = build_fold_dict(cv_distribution)

            #compute crucial and non crucial for drugs and proteins
            crucial_prots, non_crucial_prots = cc_helper(cv_distribution, protfold_dd, prots)
            crucial_drugs, non_crucial_drugs = cc_helper(cv_distribution, drugfold_dd, drugs)

            return crucial_prots, non_crucial_prots, crucial_drugs, non_crucial_drugs

        def manage_cruciality(cvdist, crucial_prots, non_crucial_prots, non_crucial_drugs):

            def find_prot(prot):
                for i,fold in enumerate(cvdist):
                    for dtis in fold:
                        if dtis[1] == prot:
                            return dtis,i

                print("ERROR: NOT FOUND")
                raise Exception

            def double_dti(foldi1, dti1, dti2, cvdist):

                #get indexes
                ind1 = cvdist[foldi1].index(dti1)

                #double
                cvdist[foldi1][ind1] = (dti1[0], dti2[1], 0)

                return cvdist

            #go through crucial_prots
            for crucial_prot in crucial_prots:

                #init temp var
                swap_done = False

                #find the fold is in
                crucial_dti, fold_prot_i = find_prot(crucial_prot)

                #go through the remaining folds
                rem_folds = set(range(foldnum)).difference(set([fold_prot_i]))

                for foldi in rem_folds:
                    
                    #set fold
                    fold = cvdist[foldi]

                    #find a DTI with no_crucial prot and no_crucial drug + negative label + drug are not the same so we can assure new label = 0
                    for dti in fold:
                        drug, prot, label = dti[0], dti[1], dti[2]

                        if not label and drug != crucial_dti[0]:

                            if drug in non_crucial_drugs and prot in non_crucial_prots:

                                cvdist = double_dti(foldi1 = foldi, 
                                                    dti1 = dti, dti2 = crucial_dti, 
                                                    cvdist = cvdist)
                                #print(f"Crucial prot {crucial_prot} doubled")
                                swap_done = True
                                break

                    if swap_done:
                        break
                                                        
            return cvdist

        while True:

            #compute pos_neg_interactions
            pos_neg_interactions = list(set(f.reduce(lambda a,b: a+b, cv_distribution)))

            #compute cruciality
            crucial_prots, non_crucial_prots, crucial_drugs, non_crucial_drugs = compute_cruciality(cv_distribution, pos_neg_interactions)

            #print(f"Crucial prots {len(crucial_prots)} - Crucial drugs {len(crucial_drugs)}")

            if not len(crucial_prots):
                #print(f"There are no crucial DTIs!\n")
                return cv_distribution
            else:
                #print(f"There are some crucial DTIs ({len(crucial_prots)})!\n")
                cv_distribution = manage_cruciality(cv_distribution, crucial_prots, non_crucial_prots, non_crucial_drugs)
 
    def Sd_St_reorder(cv_distribution, mode = 'Sd'):

        #We are going to divide DTIs in 2 groups, crucial and non-crucial
        #For a subsampling dataset to accomplish Sp requirement:

        #Drugs with +=2 targets > Drugs with == 1 target
        #(OR)
        #Proteins with +=2 drugs > Proteins with == 1 drugs

        ## Non-Crucial DTI (Di - Tj) requirement
        # Di is at least in 3 folds (2 folds would accomplish Sp requirement ✓, so if its in 3 then if moved will still accomplish)
        # Tj also accomplish this requirement

        ## Crucial DTI (Di - Tj) requirement
        # Tj is only present in 1 fold, so it needs to be swap with a Non-Crucial DTI

        def compute_cruciality(cv_distribution, pos_neg_interactions, mode = 'Sd'):

            def build_fold_dict(cv_distribution):

                def bfd_helper(mode='drug'):

                    if mode=='drug':
                        ind = 0
                    elif mode == 'prot':
                        ind = 1

                    #init default dict  (default int is 0)
                    fold_dd = dd(lambda: dd(int))

                    for i,fold in enumerate(cv_distribution):
                        fold_elements = set(list(map(lambda y: y[ind], fold)))
                        for element in fold_elements:
                            fold_dd[i][element] += 1
                        
                    return fold_dd

                protfold_dd = bfd_helper(mode='prot')
                drugfold_dd = bfd_helper(mode='drug')

                return drugfold_dd,protfold_dd

            def cc_helper(cv_distribution, fold_dd, elements):

                #init crucial and non-crucial lists
                crucial = []
                non_crucial = []

                for elm in elements:

                    prot_trues = sum([1 if fold_dd[i][elm] >= 1 else 0 for i,_ in enumerate(cv_distribution)])

                    if prot_trues >= 3:
                        non_crucial.append(elm)
                    elif prot_trues == 1:
                        crucial.append(elm)

                    # debugging
                    # if prot in [336, 682, 775, 1502]:
                    #     print(prot_trues)

                return crucial, non_crucial

            #get the proteins, as the drugs are well distributed by the way we have generated cv_distribution
            prots = set(list(map(lambda x: x[1], pos_neg_interactions)))
            drugs = set(list(map(lambda x: x[0], pos_neg_interactions)))

            #retrieve prot fold dd
            drugfold_dd, protfold_dd = build_fold_dict(cv_distribution)

            #compute crucial and non crucial for drugs and proteins
            
            crucial_prots, non_crucial_prots = cc_helper(cv_distribution, protfold_dd, prots)
            crucial_drugs, non_crucial_drugs = cc_helper(cv_distribution, drugfold_dd, drugs)

            if mode == 'Sd':
                return crucial_prots, non_crucial_prots, crucial_drugs, non_crucial_drugs
            elif mode== 'St':
                return crucial_drugs, non_crucial_drugs, crucial_prots, non_crucial_prots, 

        def manage_cruciality(cvdist, crucial_prots, non_crucial_prots, crucial_drugs, mode= 'Sd'):

            def find_crucial(elem, mode = 1):
                for i,fold in enumerate(cvdist):
                    for dtis in fold:
                        if dtis[mode] == elem:
                            return dtis,i

                print("ERROR: NOT FOUND")
                raise Exception

            def double_dti(foldi1, dti1, dti2, cvdist, mode = 'Sd'):

                #get indexes
                ind1 = cvdist[foldi1].index(dti1)

                if mode == 'Sd':
                    #double
                    cvdist[foldi1][ind1] = (dti1[0], dti2[1], 0)

                elif mode == 'St':
                    #double
                    cvdist[foldi1][ind1] = (dti2[0], dti1[1], 0)

                return cvdist

            #go through crucial_prots
            for crucial_prot in crucial_prots:

                #init temp var
                swap_done = False

                if mode == 'Sd':
                    #find the fold the crucial element is in
                    crucial_dti, fold_crucial_i = find_crucial(crucial_prot)

                elif mode == 'St':
                    #find the fold the crucial element is in
                    crucial_dti, fold_crucial_i = find_crucial(crucial_prot,0)

                #go through the remaining folds
                rem_folds = set(range(foldnum)).difference(set([fold_crucial_i]))

                for foldi in rem_folds:
                    
                    #set fold
                    fold = cvdist[foldi]

                    #find a DTI with no_crucial prot and no_crucial drug + negative label 
                    # + drug are not the same so we can assure new label = 0
                    for dti in fold:
                        drug, prot, label = dti[0], dti[1], dti[2]

                        if mode == 'Sd':

                            if not label and drug != crucial_dti[0]:

                                if drug in crucial_drugs and prot in non_crucial_prots:

                                    cvdist = double_dti(foldi1 = foldi, 
                                                        dti1 = dti, dti2 = crucial_dti, 
                                                        cvdist = cvdist)
                                                        
                                    #print(f"Crucial prot {crucial_prot} doubled")
                                    swap_done = True
                                    break

                        elif mode == 'St':

                            if not label and prot != crucial_dti[1]:

                                if prot in crucial_drugs and drug in non_crucial_prots:

                                    cvdist = double_dti(foldi1 = foldi, 
                                                        dti1 = dti, dti2 = crucial_dti, 
                                                        cvdist = cvdist, mode = mode)

                                    #print(f"Crucial drug {crucial_prot} doubled")
                                    swap_done = True
                                    break

                    if swap_done:
                        break
                                                        
            return cvdist

        while True:

            #compute pos_neg_interactions
            pos_neg_interactions = list(set(f.reduce(lambda a,b: a+b, cv_distribution)))

            #compute cruciality
            crucial_prots, non_crucial_prots, crucial_drugs, _ = compute_cruciality(cv_distribution, pos_neg_interactions, mode)

            #print(f"Crucial prots {len(crucial_prots)} - Crucial drugs {len(crucial_drugs)}")

            if not len(crucial_prots):
                #print(f"There are no crucial DTIs!\n")
                return cv_distribution
            else:
                #print(f"There are some crucial DTIs ({len(crucial_prots)})!\n")
                cv_distribution = manage_cruciality(cv_distribution, crucial_prots, non_crucial_prots, crucial_drugs, mode)
                
    def Kfold_from_lists(cv_distribution):
        for i in range(len(cv_distribution)):
            train_edges = set_to_matrix(f.reduce(lambda a,b : a+b, cv_distribution[:i] + cv_distribution[i+1:]))
            test_edges = set_to_matrix(cv_distribution[i])
            yield train_edges, test_edges

    #init seed cv list
    seed_cv_list = []

    #fix seed
    r.seed(0)
    if RMSD_dict_opt:
        n_seeds = 1
        print("Using only 1 seed for RMSD option!")
    seeds = [r.randint(1,10000) for _ in range(n_seeds)]

    print('Performing 10-CV fold for each seed')
    for seed in seeds:
        print(f"seed {seed}")
        drug_interactions_dd = get_interactions_dict(DTIs, seed, subsampling=subsampling)

        # append all interactions
        pos_neg_interactions = list(f.reduce(lambda a,b: a+b, drug_interactions_dd.values()))
        
        #check % of positives/negatives
        pos_percentage = sum(list(map(lambda x: x[2],pos_neg_interactions)))/len(pos_neg_interactions)
        print(f"Positives -> {round(pos_percentage,2)*100}%, Negatives -> {round(1-pos_percentage,2)*100} %")

        #init list to distribute edges in a Sp way.
        cv_distribution = [[] for _ in range(foldnum)]

        if mode == 'Sp':

            #both targets and drugs have been already seen during the training, but this exact DTI is new.
            for i,interaction in enumerate(pos_neg_interactions):
                #get positives for that drug
                cv_distribution[i%foldnum].append(interaction)

            #print(f"Fold sizes {list(map(len,cv_distribution))}")

            if subsampling and not RMSD_dict_opt:
                cv_distribution = Sp_reorder(cv_distribution)

        elif mode == 'Sd':

            #drugs are new, targets have been seen during the training
            pos_neg_interactions_dd = dd(list)
            for interaction in pos_neg_interactions:
                pos_neg_interactions_dd[interaction[0]].append(interaction)

            cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions_dd)

            if subsampling and not RMSD_dict_opt:
                cv_distribution = Sd_St_reorder(cv_distribution, mode= 'Sd')

        elif mode == 'St':

            pos_neg_interactions_dd = dd(list)
            #prots are new, drugs have been seen during the training
            for interaction in pos_neg_interactions:
                pos_neg_interactions_dd[interaction[1]].append(interaction)

            cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions_dd)

            if subsampling and not RMSD_dict_opt:
                cv_distribution = Sd_St_reorder(cv_distribution, mode = 'St')


        if only_distribution:
            print("Only CV distribution has been generated!")
            return cv_distribution, Drug_inv_dd, Prot_inv_dd

        #generate the interaction matrix
        pos_neg_matrix = set_to_matrix(f.reduce(lambda a,b: a+b, cv_distribution))

        ## MIKEL OPTION
        if not cvopt:

            cv_list = []

            test_df = cv_distribution[0]
            df_train, df_val = train_test_split(list(itertools.chain.from_iterable(cv_distribution[1:])), 
                                                test_size = validation_ratio/(train_ratio + validation_ratio),
                                                random_state = None, shuffle = False )

            total_len = len(df_train) + len(df_val) + len(test_df)
            
            print(f"Initial split dimensions {train_ratio}, {validation_ratio}, {test_ratio}")
            print(f"Final split dimensions {len(df_train)/total_len}, {len(df_val)/total_len}, {len(test_df)/total_len}")

            if ttv_names:

                names_train, names_val = names_from_edges(df_train, df_val)
                names_train, names_test = names_from_edges(df_train, test_df)

                cv_list = [names_train, names_val, names_test]
            
            else:

                cv_list = [df_train , df_val, test_df]


            seed_cv_list.append(cv_list)

        else:

            #init the cv list
            cv_list = []

            for train_edges, test_edges in Kfold_from_lists(cv_distribution):

                if not cvopt_no_names:

                    #--positives--
                    train_edges_pos, test_edges_pos = train_edges[train_edges[:,2] == 1,:-1], test_edges[test_edges[:,2] == 1,:-1]
                    
                    ##create names matrix from edges list
                    names_train_pos, names_test_pos = names_from_edges(train_edges_pos, test_edges_pos)

                    #--negatives--
                    train_edges_neg, test_edges_neg = train_edges[train_edges[:,2] == 0,:-1], test_edges[test_edges[:,2] == 0,:-1]
                    
                    ##create names matrix from edges list
                    names_train_neg, names_test_neg = names_from_edges(train_edges_neg,test_edges_neg)

                    #print(f"Train pos {len(names_train_pos)}, Train neg {len(names_train_neg)}, Test pos {len(names_test_pos)}, Test neg {len(names_test_neg)}")

                    #add each fold
                    cv_list.append((names_train_pos,names_train_neg,names_test_pos,names_test_neg))

                else:

                    cv_list.append((train_edges,test_edges))

            #add each group of folds for each seed
            seed_cv_list.append(cv_list)


    if RMSD_dict_opt:

        #Get negative proteins
        neg_prots = list(set(map(lambda x: x[1], filter(lambda x: x[2] == 0, pos_neg_interactions))))

        #Get positive interactions from the pos_neg_interactions and positive proteins for the final fold, make sure they are not new
        pos_prots = list(set(map(lambda x: Prot_inv_dd[x[1]], filter(lambda x: x[2] == 1, pos_neg_interactions))))
        pos_prots_final = list(set(map(lambda x: x[1], filter(lambda x: x[2] == 1, positive_final_fold))))

        assert not set(pos_prots_final).difference(set(pos_prots))

        negative_final_fold = r.sample(negative_final_fold, len(positive_final_fold))
        final_fold = positive_final_fold + negative_final_fold
        prot_info_dict = {'neg_prot_dict': Counter(neg_prots),
                          'neg_percentage': round(len(neg_prots) / Prot_L * 100,2),
                          'final_fold' : final_fold}
        print(f"{len(positive_final_fold)} positives and {len(negative_final_fold)} negatives selected for final fold, negative is using a {prot_info_dict['neg_percentage']}% of total proteins")

        return seed_cv_list, prot_info_dict

    else:

        return seed_cv_list

#Lets use Yamanishi NR as an example
##Load dataset
#fpath = os.path.join(os.getcwd(),'DB','Data','Yamanashi_et_al_GoldStandard','IC','interactions','ic_admat_dgc_mat_2_line.txt')
#fpath = os.path.join(os.getcwd(),'DB','Data','Yamanashi_et_al_GoldStandard','NR','interactions','nr_admat_dgc_mat_2_line.txt')
#fpath = os.path.join(os.getcwd(),'DB','Data','Davis_et_al','tdc_package_preprocessing','DAVIS_et_al_2line.tsv')
#fpath = os.path.join(os.getcwd(), 'DB', 'Data', 'BIOSNAP', 'ChG-Miner_miner-chem-gene', 'ChG-Miner_miner-chem-gene.tsv')
#fpath = os.path.join(os.getcwd(), 'DTINet', 'Data', 'BindingDB', 'final_dtis_BindingDB.tsv')
#DTIs = pd.read_csv(fpath, sep='\t') ## MAKE SURE THE HEADER OPTION IS ON/OFF DEPENDING ON THE DATASET!
#DTIs.columns = ['Drug', 'Protein']

#check splits
def check_splits(splits, verbose=False, foldnum=10):

    def check_condition():
        
        for seed in range(len(splits)):
            print(f"Seed {seed}")
            drugs_counter, prots_counter = 0,0
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
                if drugs_counter == (foldnum):
                    print("Sd split configuration accomplished!")

                else:
                    print(f"Sd split configuration not exactly accomplished! ({drugs_counter})")

            if prots_counter > 0:
                if prots_counter == (foldnum):
                    print("St split configuration accomplished!")
                else:
                    print(f"St split configuration not exactly accomplished! ({prots_counter})")

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


# ------------------------------------------------------- Sp ------------------------------------------------------------------ #
##Get 5-seed 10-fold CV Sp (all nodes are seeing during the training)
# sp_splits, prot_info_dict = generate_splits(DTIs, mode= 'Sp', subsampling=True, foldnum=10, 
#                             negative_to_positive_ratio = 1, cvopt=True, RMSD_threshold = 10,
#                             RMSD_dict_opt=True, include_diagonal_RMSD=False)

#FOR GUILLE - LINE 86 (SP)
#cv_distr, inv_drug_dd, inv_prot_dd = generate_splits(DTIs, mode= 'Sp', subsampling=False, foldnum=10, RMSD_dict_opt=False, only_distribution=True)

##check sp split
#check_splits(sp_splits, verbose=False) 

##check distribution
#print_cv_distribution(DTIs, sp_splits[0])

# ------------------------------------------------------- Sd ------------------------------------------------------------------- #
##Get 5-seed 10-fold CV Sd (some drugs are not seeing during the training)
#sd_splits = generate_splits(DTIs, mode= 'Sd', subsampling=True, foldnum=10, RMSD_dict_opt=True)

##check splits
#check_splits(sd_splits,verbose=False) 

##check distributions
#print_cv_distribution(DTIs,sd_splits)
# ------------------------------------------------------- St ------------------------------------------------------------------- #
##Get 5-seed 10-fold CV Sd (some targets are not seeing during the training)
#st_splits = generate_splits(DTIs, mode= 'St', subsampling=True, foldnum=10, RMSD_dict_opt=True)

##check splits
#check_splits(st_splits,verbose=False) 