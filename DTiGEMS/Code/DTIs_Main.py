# coding: utf-8
# All needed packages
import argparse
import node2vec
from gensim.models import Word2Vec
import pandas as pd
import math as math
import numpy as np
import csv
import time
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import  RandomForestClassifier, AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from imblearn.over_sampling import RandomOverSampler, SMOTE, ADASYN
from sklearn.metrics import *
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import MaxAbsScaler

# Import my files
from n2v_mainFunctions import *
from load_datasets import *
from pathScores_functions import *
from get_Embeddings_FV import *
from training_functions import *
from GIP import *
from snf_code import *
######################################## START MAIN #########################################
#############################################################################################

def main():
    # get the parameters from the user
    args = parse_args()
    ## get the start time to report the running time
    t1 = time.time()

    ### Load the input data - return all pairs(X) and its labels (Y)..
    allD, allT, allDsim, allTsim, DrTr, R, X, Y = load_datasets()

    # create 2 dictionaries for drugs. First one the keys are their order numbers
    #the second  one the keys are their names -- same for targets
    drugID = dict([(d, i) for i, d in enumerate(allD)])
    targetID = dict([(t, i) for i, t in enumerate(allT)])
    #-----------------------------------------
    ###### Define different classifiers

    #Random Forest...........
    rf = RandomForestClassifier(n_estimators=200 ,n_jobs=10,random_state= 55, class_weight='balanced', criterion='gini', max_features = 'auto') 

    # Adaboost .......
    ab = AdaBoostClassifier(DecisionTreeClassifier(criterion='gini', splitter='best', max_depth=15, min_samples_split=2,
                        min_samples_leaf=1, max_features=None, random_state=1,max_leaf_nodes=None, 
                        class_weight= 'balanced' ), algorithm="SAMME", n_estimators=100,random_state=32)

    # 10-folds Cross Validation...............
    skf = StratifiedKFold(n_splits=10, shuffle = True, random_state = 22)
    skf.get_n_splits(X, Y)
    # all evaluation listsSelect..
    correct_classified = []
    ps = []
    recall = []
    roc_auc = []
    average_precision = []
    f1 = []
    Pre = []
    Rec = []
    AUPR_TEST = []
    TN = []
    FP = []
    FN = []
    TP = []

    foldCounter = 1     # fold counter
    for train_index, test_index in  skf.split(X,Y):
        
        print("*** Working with Fold %i :***" %foldCounter)

        # create known interaction of training part 
        R_train_pos = []
        for i in train_index:
            if(Y[i]==1):
                dr = X[i,0]
                tr = X[i,1]
                Rpos_data = dr, tr, Y[i]
                R_train_pos.append(Rpos_data)

        R_tr = pd.DataFrame(R_train_pos)

        #first thing with R train to remove all edges in test (use it when finding path)
        train_DT_Matrix = Mask_test_index(test_index, X, Y, DrTr, drugID, targetID)
        DrTr_train = train_DT_Matrix.transpose()

        #^^^^^^^^^^^^^^^^^^^^^^^^^ GIP SIMILARITY ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ## get GIP Similarity from training known interactions, TT sim, DD sim
        DT_impute_D = impute_zeros(DrTr_train,allDsim[1])
        DT_impute_T = impute_zeros(np.transpose(DrTr_train),allTsim[1])

        GIP_D = Get_GIP_profile(np.transpose(DT_impute_D),"d")
        GIP_T = Get_GIP_profile(DT_impute_T,"t")

        # ^^^^^^^^^^^^^^^^^^^^^ SNF ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #fused similarity matrices into one matrix using SNF

        DDsim = []
        TTsim = []

        #DDsim.append(normDD_Sim)
        DDsim.append(GIP_D)
        for sim in allDsim:
                DDsim.append(sim)
        
        #TTsim.append(normTT_Sim)
        TTsim.append(GIP_T)      
        for sim in allTsim:
                TTsim.append(sim)

        fused_simDr = SNF(DDsim,K=5,t=2,alpha=1.0)
        fused_simTr = SNF(TTsim,K=5,t=2,alpha=1.0)

        #^^^^^^^^^^^^^^^^^^^^^^^^^DTIs with node2vec code^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        k = 5
        # remove all weak edges by keep strongest k similarities
        TSsim = Strongest_k_sim(allTsim[0] ,k)
        DSsim = Strongest_k_sim(allDsim[0],k)

        # convert GIP similarity matrix into edgelist
        PP_edgeList = edgeList(TSsim , allT)
        DD_edgeList = edgeList(DSsim, allD )
        # create the edge list of the whole graph including (DD, PP, R_train)
        frames = [R_tr, DD_edgeList, PP_edgeList]
        allEdgeList = pd.concat(frames)
        
        # write it in a file for next reading of node2vec
        input_filename = 'edgeLists/edgeList_%i.txt '%foldCounter
        allEdgeList.to_csv(input_filename, header = None, index = None, sep=' ')

        # write edgeList and then read it for NODE2VEC generate Embedding
        # read the graph to feed it into node2vec
        train_graph_file = open(input_filename, 'r')
        args.input = train_graph_file  
        print('input is',args.input)

        #<<<<<<<<<<<<<<<<<<<<<<<<><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # apply n2v on R positive training part and get the embeddings (protien_drugs FV)
        nx_G = read_graph(args)
        G = node2vec.Graph(nx_G, args.directed, args.p, args.q)
        G.preprocess_transition_probs()
        walks = G.simulate_walks(args.num_walks, args.walk_length)
        learn_embeddings(walks, args)

        # get the feature vector of R_positive training part (think about sending part!!??)
        targetFV, drugFV = get_FV_drug_target(args, foldCounter, allT, allD)
        
        # Calculate cosine similarity for each drug pair, and for each target pair
        cos_simDD = Cosine_Similarity(drugFV)
        cos_simTT = Cosine_Similarity(targetFV)
        # normalize simiarities to be in positive range [0,1]
        cos_simDD = normalizedMatrix(cos_simDD)
        cos_simTT  = normalizedMatrix(cos_simTT )
        #---------------------------------------------------------------------   

        # Generate all featres from the matrix multiplication of each path strucutre
        # list for each feature (Graph G1)
        sumDDD, maxDDD = DDD_TTT_sim(fused_simDr)
        sumTTT, maxTTT= DDD_TTT_sim(fused_simTr)
        
        sumDDT,maxDDT = metaPath_Dsim_DT(fused_simDr,DrTr_train,2,True) 
        sumDTT,maxDTT = metaPath_DT_Tsim(fused_simTr,DrTr_train,2,True)

        sumDDDT,_ = metaPath_Dsim_DT(sumDDD,DrTr_train,3,True)
        _,maxDDDT = metaPath_Dsim_DT(maxDDD,DrTr_train,3,True)

        sumDTTT,_ = metaPath_DT_Tsim(sumTTT,DrTr_train,3,True)
        _,maxDTTT = metaPath_DT_Tsim(maxTTT,DrTr_train,3,True)

        sumDTDT,maxDTDT = metaPath_DTDT(DrTr_train)
        sumDDTT,maxDDTT = metaPath_DDTT(DrTr_train,fused_simDr,fused_simTr, True)
    #============================================================================== 
        # Generate all featres from the matrix multiplication of each path strucutre
        # list for each feature (Graph G2)
        sumDDD2, maxDDD2 = DDD_TTT_sim(cos_simDD)
        sumTTT2, maxTTT2= DDD_TTT_sim(cos_simTT)
        
        sumDDT2,maxDDT2 = metaPath_Dsim_DT(cos_simDD,DrTr_train,2,True) 
        sumDTT2,maxDTT2 = metaPath_DT_Tsim(cos_simTT,DrTr_train,2,True)

        sumDDDT2,_ = metaPath_Dsim_DT(sumDDD2,DrTr_train,3,True)
        _,maxDDDT2 = metaPath_Dsim_DT(maxDDD2,DrTr_train,3,True)

        sumDTTT2,_ = metaPath_DT_Tsim(sumTTT2,DrTr_train,3,True)
        _,maxDTTT2 = metaPath_DT_Tsim(maxTTT2,DrTr_train,3,True)

        sumDTDT2,maxDTDT2 = metaPath_DTDT(DrTr_train)
        sumDDTT2,maxDDTT2 = metaPath_DDTT(DrTr_train,cos_simDD,cos_simTT, True)
        
    #============================================================================== 
    ### Build feature vector and class labels
        DT_score = []
        for i in range(len(allD)):
            for j in range(len(allT)):        
                pair_scores = (allD[i], allT[j],\
                            # path scores from G1
                               sumDDT[i][j],sumDDDT[i][j],\
                               sumDTT[i][j],sumDTTT[i][j], sumDDTT[i][j], sumDTDT[i][j],\
                               maxDDT[i][j],maxDDDT[i][j], \
                               maxDTT[i][j],maxDTTT[i][j],maxDDTT[i][j],maxDTDT[i][j],\
                            # path scores from G2
                               sumDDT2[i][j],sumDDDT2[i][j],\
                               sumDTT2[i][j],sumDTTT2[i][j], sumDDTT2[i][j], sumDTDT2[i][j],\
                               maxDDT2[i][j],maxDDDT2[i][j], \
                               maxDTT2[i][j],maxDTTT2[i][j],maxDDTT2[i][j],maxDTDT2[i][j])
                DT_score.append(pair_scores)
        
        features = []
        class_labels = []
        DT_pair = []
        # Build the feature vector - Concatenate features from G1,G2
        for i in range(len(DT_score)):
            dr = DT_score[i][0]
            tr = DT_score[i][1] 
            edgeScore = DT_score[i][2], DT_score[i][3], DT_score[i][4],DT_score[i][5],DT_score[i][6],\
            DT_score[i][8],DT_score[i][9], DT_score[i][10], DT_score[i][11],DT_score[i][12],\
            DT_score[i][14], DT_score[i][15],DT_score[i][16],DT_score[i][17],DT_score[i][18],\
            DT_score[i][20],DT_score[i][21],DT_score[i][22],DT_score[i][23]#,DT_score[i][24]

            dt = DT_score[i][0], DT_score[i][1]
            DT_pair.append(dt)
            features.append(edgeScore)
            # same label as the begining
            label = R[dr][tr]
            class_labels.append(label)

        ## Start Classification Task
        # featureVector and labels for each pair
        XX = np.asarray(features)
        YY = np.array(class_labels)

        #fit the model
        max_abs_scaler = MaxAbsScaler()
        X_train = max_abs_scaler.fit(XX[train_index]) 
        X_train_transform = max_abs_scaler.transform(XX[train_index])

        X_test_transform = max_abs_scaler.transform(XX[test_index])

        rf.fit(X_train_transform, YY[train_index])
        predictedClass= rf.predict(X_test_transform)
        predictedScore = rf.predict_proba(X_test_transform)[:, 1]

        print("@@ Validation and evaluation of fold %i @@" %foldCounter)

        #print('Done.\nAccuracy: %f' % accuracy_score(Y[test_index], predictedClass))
        cm = confusion_matrix(YY[test_index], predictedClass)
        TN.append(cm[0][0])
        FP.append(cm[0][1])
        FN.append(cm[1][0])
        TP.append(cm[1][1])
        print("Confusion Matrix for this fold")
        print(cm)

        print("Correctly Classified Instances: %d" %accuracy_score(Y[test_index], predictedClass, normalize=False))
        correct_classified.append(accuracy_score(Y[test_index], predictedClass, normalize=False))

        #print("Precision Score: %f" %precision_score(Y[test_index], predictedClass))
        ps.append(precision_score(Y[test_index], predictedClass,average='weighted'))

        #print("Recall Score: %f" %recall_score(Y[test_index], predictedClass)
        recall.append(recall_score(Y[test_index], predictedClass, average='weighted'))

        print("Area ROC: %f" %roc_auc_score(Y[test_index], predictedClass))
        roc_auc.append(roc_auc_score(Y[test_index], predictedScore))

        print("F1 Score: %f" %f1_score(Y[test_index], predictedClass, average='weighted'))
        f1.append(f1_score(Y[test_index], predictedClass))

        p, r, _ = precision_recall_curve(Y[test_index],predictedScore,pos_label=1)
        aupr = auc(r, p)
        AUPR_TEST.append(aupr)

        print("AUPR auc(r,p) = %f" %aupr)
        #print("AUPR = %f" %average_precision_score(Y[test_index],  predictedClass))
        average_precision.append(average_precision_score(Y[test_index], predictedClass))

        print(classification_report(Y[test_index], predictedClass))
        print('-------------------------------------------------------------')
        foldCounter += 1
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

    ############# Evaluation Metrics ####################################
    # Confusion matrix for all folds
    ConfMx = np.zeros((cm.shape[0],cm.shape[0]))
    ConfMx[0][0] = str( np.array(TN).sum() )
    ConfMx[0][1] = str( np.array(FP).sum() )
    ConfMx[1][0] = str( np.array(FN).sum() )
    ConfMx[1][1] = str( np.array(TP).sum() )

    ### Print Evaluation Metrics using all and each path category scores FV,  RF classifier, and NR dataset..
    print("Result(Correct_classified): " + str( np.array(correct_classified).sum() ))
    print("Results:precision_score = " + str( np.array(ps).mean().round(decimals=3) ))
    print("Results:recall_score = " + str( np.array(recall).mean().round(decimals=3) ))
    print("Results:roc_auc = " + str( np.array(roc_auc).mean().round(decimals=3) ))
    print("Results:f1 = " + str( np.array(f1).mean().round(decimals=3) ))
    print("Results: AUPR on Testing auc(r,p) = " + str( np.array(AUPR_TEST).mean().round(decimals=3)))
    print("Confusion matrix for all folds")
    print(ConfMx)
    print('_____________________________________________________________')
    print('Running Time for the whole code:', time.time()-t1)  
    print('_____________________________________________________________')  
#####-------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
#####---------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ######################################################################