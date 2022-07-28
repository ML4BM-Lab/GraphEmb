import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from collections import defaultdict as dd
import random as r
import functools as f
from tqdm import tqdm
import itertools
from tqdm.contrib.itertools import product
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from tqdm import tqdm
import logging
import os
import helper_functions_dtinet as hf
import argparse
from tqdm.contrib.itertools import product
from sklearn.model_selection import train_test_split



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

	def get_interactions_dict(DTIs, seed, subsampling, RMSD_dict = None):

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

		#init default dict (list)
		interactions_dd = dd(list)

		#get prng
		prng = np.random.RandomState(seed)

		#get positives
		for d,p in zip(DTIs['Drug'],DTIs['Protein']):
			interactions_dd[Drug_dd[d]].append((Drug_dd[d],Prot_dd[p],1))
			   
		#add negatives (subsample to have 50%-50%)
		#go through all drugs/proteins
		for i, elementid in enumerate(interactions_dd):
			
			#print(f"elementid {elementid} in i {i}")
			#subsample from the negatives and add it to interactions dictionary
			#drugs if swap = False | proteins if swap = True
			pos_element = list(map(lambda x: x[1], interactions_dd[elementid])) #x[1] = prot
			neg_element = set(range(Prot_L)).difference(set(pos_element))
			#print(f"Positive element {len(pos_element)}, negative elements {len(neg_element)}")

			if subsampling:
				#check if we can subsample all
				if RMSD_dict is None:
					logging.debug('rmsd none in script')
					neg_sampled_element = r.sample(neg_element,min(len(neg_element),len(pos_element))) #50%-50% (modify if different proportions are desired)
				else:
					logging.debug('rmsd selected in script')
					neg_sampled_element = get_targets_for_drugs_RMSD(pos_element, neg_element, RMSD_dict, Prot_inv_dd)
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

		#generate names matrix from positive edges
		train_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p]) for d,p in train_edges]

		#same with negatives
		test_names_matrix = [(Drug_inv_dd[d], Prot_inv_dd[p]) for d,p in test_edges]

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

	for seed in tqdm([7183, 556, 2, 81, 145], desc='Performing 10-CV fold for each seed'):

		drug_interactions_dd = get_interactions_dict(DTIs, seed, subsampling=subsampling, RMSD_dict=RMSD_dict)

		# append all interactions
		pos_neg_interactions = list(f.reduce(lambda a,b: a+b, drug_interactions_dd.values()))
		
		#check % of positives/negatives
		pos_percentage = sum(list(map(lambda x: x[2],pos_neg_interactions)))/len(pos_neg_interactions)
		#print(f"Positives -> {round(pos_percentage,2)*100}%, Negatives -> {round(1-pos_percentage,2)*100} %")

		#init list to distribute edges in a Sp way.
		cv_distribution = [[] for _ in range(foldnum)]

		if mode == 'Sp':

			#both targets and drugs have been already seen during the training, but this exact DTI is new.
			for i,interaction in enumerate(pos_neg_interactions):
				#get positives for that drug
				cv_distribution[i%foldnum].append(interaction)

			#print(f"Fold sizes {list(map(len,cv_distribution))}")

			if subsampling:
				cv_distribution = Sp_reorder(cv_distribution)

		elif mode == 'Sd':

			#drugs are new, targets have been seen during the training
			pos_neg_interactions_dd = dd(list)
			for interaction in pos_neg_interactions:
				pos_neg_interactions_dd[interaction[0]].append(interaction)

			cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions_dd)

			if subsampling:
				cv_distribution = Sd_St_reorder(cv_distribution, mode= 'Sd')

		elif mode == 'St':

			pos_neg_interactions_dd = dd(list)
			#prots are new, drugs have been seen during the training
			for interaction in pos_neg_interactions:
				pos_neg_interactions_dd[interaction[1]].append(interaction)

			cv_distribution = optimize_folds(cv_distribution, pos_neg_interactions_dd)

			if subsampling:
				cv_distribution = Sd_St_reorder(cv_distribution, mode = 'St')

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

			#init the cv list
			cv_list = []

			for train_edges, test_edges in Kfold_from_lists(cv_distribution):
			#for train_edges, test_edges in Kfold_from_lists([[1],[2],[3],[4],[5]]):
				# print(train_edges)
				# print(test_edges)
				# print()
				  
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

	for i,gen in enumerate(tqdm(names, desc='rmsd dict')):

		#get genes and names
		rmsd_i = RMSD[i,0:i].tolist() + RMSD[i,i+1:].tolist()
		names_i = names[0:i] + names[i+1:]

		#add the entry
		RMSDdict[gen] = dict(zip(names_i,rmsd_i))


	return RMSDdict

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





def get_idx_matlab(wdir, sp_splits):
	# for dtinet we need the tuple as (protein, drug)
	# make translation as function --> *
	index_protein_file = 'protein.txt'
	index_drug_file = 'drug.txt'
	# def translate_index():
	drugs_list = np.loadtxt(os.path.join(wdir, index_drug_file), dtype=str).tolist()
	proteins_list = np.loadtxt(os.path.join(wdir, index_protein_file), dtype=str).tolist()
	# build dics
	drug_pos_d = dict(zip(drugs_list,list(range(len(drugs_list)))))
	protein_pos_d = dict(zip(proteins_list,list(range(len(proteins_list)))))
	# loop over seeds, folds, set_type 
	# shows total of the tqdm object 
	nseed, nfold, set_type = 0,0,0 # init
	for nseed, nfold, set_type in product(range(len(sp_splits)), range(len(sp_splits[nseed])), range(len(sp_splits[nseed][nfold])),
										desc='changing ID for matrix index'): 
		# iteration over pairs to change tuples in list
		for pair in range(len(sp_splits[nseed][nfold][set_type])): # translate pairs
			# this is  drug_protein
			drug_coo = drug_pos_d[sp_splits[nseed][nfold][set_type][pair][0]]  
			prot_coo = protein_pos_d[sp_splits[nseed][nfold][set_type][pair][1]] 
			idx_matlab = np.ravel_multi_index((drug_coo,prot_coo), (len(drugs_list), len(proteins_list)), order='F') + 1 # matlab +1
			# add here to change between drug index target index to matlab intex
			sp_splits[nseed][nfold][set_type][pair] = (idx_matlab) # uncoment
			#sp_splits[nseed][nfold][set_type][pair] = (drug_coo+1, prot_coo+1)
	return sp_splits



def main():
	'''
	generate one to one index foldere & files
	considering splits and subsampling 
	'''
	parser = argparse.ArgumentParser() 
	parser.add_argument("-v", "--verbose", dest="verbosity", action="count", default=4,
					help="Verbosity (between 1-4 occurrences with more leading to more "
						"verbose logging). CRITICAL=0, ERROR=1, WARN=2, INFO=3, "
						"DEBUG=4")
	parser.add_argument("-dbPath","--dbPath", help="Path to the database output ('BIOSNAP', 'BindingDB', 'Davis_et_al', 'DrugBank_FDA', 'E', 'GPCR', 'IC', 'NR')", type=str)
	parser.add_argument("-split_type", "--split_type", help="Select the type of split ['Sp', 'Sd', 'St'] to generate oneTooneIndex folder", type=str)
	parser.add_argument("-subsampling", help="Flag for subsampling True", action='store_true')
	parser.add_argument("-rmsd", help="Flag for RSMD True", action='store_true')

	args = parser.parse_args()

	# -) log info; 
	# define logging level
	log_levels = {
		0: logging.CRITICAL,
		1: logging.ERROR,
		2: logging.WARN,
		3: logging.INFO,
		4: logging.DEBUG,
	}
	# set the logging info
	level= log_levels[args.verbosity]
	fmt = '[%(levelname)s] %(message)s'
	logging.basicConfig(format=fmt, level=level)

	#######
	### log output detals
	logging.info(
		'''
		This script generates splits
		'''
		)
	# OUTPUT DIRECTORY
	# sanity check
	DB_PATH = args.dbPath
	SPLIT_TYPE = args.split_type
	if SPLIT_TYPE not in  ['Sp', 'Sd', 'St']:
		raise NameError('Need to specify a valid split type')

	SUBSAMPLING_TYPE = args.subsampling

	RMSD_SET = args.rmsd

	logging.info(f'Working in output folder for: {DB_PATH}')
	logging.info(f'Creating Split: {SPLIT_TYPE}')
	logging.info(f'Subsampling: {SUBSAMPLING_TYPE}')
	logging.info(f'RMSD: {RMSD_SET}')

	db_name = hf.get_DB_name(DB_PATH)
	hf.check_and_create_folder(db_name)
	wdir = os.path.join('../Data', db_name)

	##########################################

	# create index folder if it does not exists
	# select mode
	if SUBSAMPLING_TYPE:
		sub = 'subsampling'
	else:
		sub = 'nosubsampling'

	# generate RMSD dictionary
	if RMSD_SET:
		RMSD_dict = genRMSDdict()
		rmsd_name = '_rmsd'
		logging.debug('Set rmsd')
	else:
		RMSD_dict = None
		rmsd_name = ''
		logging.debug('rmsd not set')


	# create folders
	path_folder = os.path.join(wdir, f'Index_{SPLIT_TYPE}_{sub}{rmsd_name}')
	if not os.path.exists(path_folder):
		os.makedirs(path_folder)
	
	# LOAD DTIs
	fpath = os.path.join(os.getcwd(), wdir, f'final_dtis_{DB_PATH}.tsv')
	DTIs = pd.read_csv(fpath,sep='\t',header=None)
	DTIs.columns = ['Drug','Protein']
	#logging.debug(DTIs.head())


	## GENERATE SPLITS
	splits = generate_splits(DTIs, mode= SPLIT_TYPE, subsampling=SUBSAMPLING_TYPE, foldnum=10, RMSD_dict=RMSD_dict)
	
	#logging.debug(splits[0][0][0][0])
	# Convert to Matlab index type
	# this also changes (drug, protein) to (protein, drug)
	splits_matlab = get_idx_matlab(wdir, splits)
	#logging.debug(splits_matlab[0][0][0][0])
	# Save splits as .txt
	nseed, nfold = 0, 0 
	for nseed, nfold in product(range(len(splits_matlab)), range(len(splits_matlab[nseed]))):
		np.savetxt(os.path.join(path_folder, f'train_pos_{nseed+1}_{nfold+1}.txt'), splits_matlab[nseed][nfold][0], fmt='%i', delimiter=" ")
		np.savetxt(os.path.join(path_folder, f'train_neg_{nseed+1}_{nfold+1}.txt'), splits_matlab[nseed][nfold][1], fmt='%i', delimiter=" ")
		np.savetxt(os.path.join(path_folder, f'test_pos_{nseed+1}_{nfold+1}.txt'), splits_matlab[nseed][nfold][2], fmt='%i', delimiter=" ")
		np.savetxt(os.path.join(path_folder, f'test_neg_{nseed+1}_{nfold+1}.txt'), splits_matlab[nseed][nfold][3], fmt='%i', delimiter=" ")
	
#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":
	main()
#####-------------------------------------------------------------------------------------------------------------
####################### END OF THE CODE ##########################################################################
