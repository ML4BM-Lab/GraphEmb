import sys
from Graph_utils import *
from functions import cross_validation,mean_confidence_interval
from Classify import *
from SNF import *
from datetime import datetime
import numpy as np



def get_features_per_fold(R_train,D_sim,T_sim, pair):
	
	accum_DDD,max_DDD = get_two_hop_similarities(D_sim,D_sim)
	accum_TTT,max_TTT = get_two_hop_similarities(T_sim,T_sim)



	accum_DDT,max_DDT = get_drug_relation(D_sim,R_train) 
	accum_DTT,max_DTT = get_relation_target(T_sim,R_train)


	accum_DDDT,_ = get_drug_relation(accum_DDD,R_train)
	_,max_DDDT = get_drug_relation(max_DDD,R_train)


	accum_DTTT,_ = get_relation_target(accum_TTT,R_train)
	_,max_DTTT = get_relation_target(max_TTT,R_train)


	accum_DTDT,max_DTDT = get_DTDT(R_train)

	accum_DDTT,max_DDTT = get_DDTT(R_train,D_sim,T_sim)

	features = []


	
	features.append(mat2vec(accum_DDT))
	features.append(mat2vec(max_DDT))
	features.append(mat2vec(accum_DTT))
	features.append(mat2vec(max_DTT))

	features.append(mat2vec(accum_DDDT))

	features.append(mat2vec(max_DDDT))

	if pair:
		features.append(mat2vec(accum_DTDT))
		features.append(mat2vec(max_DTDT))

	
	features.append(mat2vec(accum_DTTT))

	features.append(mat2vec(max_DTTT))

	features.append(mat2vec(accum_DDTT))

	features.append(mat2vec(max_DDTT))

	return features





def get_similarities(sim_file,dMap):

	sim = []

	for line in open(sim_file).readlines():
		edge_list = get_edge_list(line.strip())
		new_matrix = make_sim_matrix(edge_list,dMap)
		#print "shape of new matrix "+str(new_matrix.shape)
		sim.append(new_matrix)
	return sim


#### main script#################
R_all_train_test = sys.argv[1] #"nr_admat_dgc_mat_2_line.txt"

D_sim_file = sys.argv[2] #"nr_D_similarities.txt"
T_sim_file = sys.argv[3] # "nr_T_similarities.txt"




(D,T,DT_signature,aAllPossiblePairs,dDs,dTs,diDs,diTs) = get_All_D_T_thier_Labels_Signatures(R_all_train_test)

R = get_edge_list(R_all_train_test)
DT = get_adj_matrix_from_relation(R,dDs,dTs)


D_sim = get_similarities(D_sim_file,dDs)
T_sim = get_similarities(T_sim_file,dTs)



row,col = DT.shape
#print "number of rows "+str(row)
#print "number of cols "+str(col)

all_matrix_index = []

for i in range(row):
	for j in range(col):
		all_matrix_index.append([i,j])


for mode in ["p","D","T"]:
	seeds = [7771, 8367, 22, 1812, 4659]


	if mode == "p":
		cv_data = cross_validation(DT,seeds,1,10)
		pair = True
	elif mode == "D":
		cv_data = cross_validation(DT,seeds,0,10)
		pair = False
	elif mode == "T":
		pair = False
		cv_data = cross_validation(np.transpose(DT),seeds,0,10)

	labels = mat2vec(DT)
	#print(labels)

	test_idx = []
	trails_AUPRs = []
	trails_AUCs = []
	trails_recall  = []
	trails_precision = []
	trails_f1 = []

	for seed in seeds:
		print "seed -> "+str(seed)+" time: "+str(datetime.now())
		total = 0
		folds_features = []
		for fold in cv_data[seed]:
			#print "a"

			if mode == "T":
				R_train = mask_matrix(DT,fold[1],True)
			else:
				R_train = mask_matrix(DT,fold[1]) #by default transpose is false
			

			DT_impute_D = impute_zeros(R_train,D_sim[0])
			DT_impute_T = impute_zeros(np.transpose(R_train),T_sim[0])
			
			GIP_D = Get_GIP_profile(np.transpose(DT_impute_D),"d")
			GIP_T = Get_GIP_profile(DT_impute_T,"t")


			WD = []
			WT = []

			for s in D_sim:
				WD.append(s)
			WD.append(GIP_D)

			for s in T_sim:
				WT.append(s)
			WT.append(GIP_T)
			D_SNF = SNF(WD,3,2)
			T_SNF = SNF(WT,3,2)

			DS_D = FindDominantSet(D_SNF,5)
			DS_T = FindDominantSet(T_SNF,5)

			np.fill_diagonal(DS_D,0)
			np.fill_diagonal(DS_T,0)
			
			
			features = get_features_per_fold(R_train,DS_D,DS_T, pair)
			folds_features.append(zip(*features))
			if mode == "T":
				test_idx.append([j*col+i for (i,j) in fold[1]])
			else:
				test_idx.append([i*col+j for (i,j) in fold[1]])

		#print total
		results = run_classification(folds_features,labels,test_idx)
		trails_AUPRs.extend(results[2])
	aupr,c1 = mean_confidence_interval(trails_AUPRs)
	print "################Results###################"
	print "Mode: %s"%mode
	print "Average AUPR: %f"%aupr
	print "std AUPR %f"%np.std(results[2])
	print "Average AUC: %f"%np.mean(results[-1])
	print "std AUC: %f"%np.std(results[-1])
	print "###########################################"
