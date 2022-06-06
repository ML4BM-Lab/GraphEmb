##We changed 3 lines in DDR, one in the demo (to print the AUC) and 
## we also change 3 other lines in the demo for inputing parameters instead of hardcoding.
##two in the SNF to be able to compute the pairwise euclidean distance

#------------------------------------ FIRST ------------------------------------------------------#
#FIRST, copy all the symmat_nonames (from PreSNF folder) folder to :/ddr/DATASETS/

# path -> /home/jfuente/data/jfuente/DTI/Input4Models/DDR/Data/{dataset}/Formatted/PreSNF
## -> docker cp symmat_nonames $dockername:/ddr/DATASETS/{dataset}

## Now run the SNF script from the folder containing all the symmat_nonames files
## -> python ../../../Similarity_selection_code/Similarity_Selection.py snflist_drugs/prots.txt (both) 0.6 0.7 (default parameters as shown in the paper)


#---------------------------------- SECOND ----------------------------------------------------#
## copy the names into the PostSNF matrix (create two files called selected_drugs.txt and selected_prots.txt)

# now copy the corresponding matrices from the 2 line folder of preSNF to the PostSNF and create two files 
# containing the SNF matrices but with the 2 line format (selected_drugs_post.txt and selected_prots_post.txt)

# copy also the DTI file from the DB folder

#------------------------------------- THIRD --------------------------------------------------#
# copy the PostSNF folder to the docker 
# docker cp PostSNF $dockername:/ddr/DATASETS/{dataset}

# run the nr_demo.py with the DTI and the drug similarity and prot similarity matrices
# nohup python ../../../Similarity_selection_code/DDR_demo/nr_demo.py {DTIs_matrix} {drug_sim} {prot_sim} > results_{dataset_name}.out &

#-------------------------------------- FOURTH ------------------------------------------------- #

# copy the result to the DDR/Results/ folder
# enjoy