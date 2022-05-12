#!/bin/bash          
# to be executed in Code folder! 

# this XXX for XXX
# change for your XXX
# copy into your machine
# docker cp /home/uveleiro/data/jfuente/DTI/Input4Models/EEG-DTI/Code/process_eegdti_matlab.sh matlab_env_t2:/DTINet


wdir=Davis_et_al # define for each database -> parse this information
folder_path=../Data/$wdir

######### sevenNets  #########
[ ! -d $folder_path/sevenNets ] && mkdir $folder_path/sevenNets
cp $folder_path/mat_*.txt $folder_path/sevenNets

### add here oneTooneIndex
## ----------------------------------------->> ************

######### sim_network #########
# data for similarity == Luo compute sim, retrieve only \network
[ ! -d $folder_path/tmp_data ] && mkdir $folder_path/tmp_data

cp $folder_path/mat_*.txt $folder_path/tmp_data
cp $folder_path/Similarity_Matrix_*.txt $folder_path/tmp_data

# now we need to copy this to a docker
# if this is not your case
# execute matlab here

# copy to docker
container=matlab_env_t2 # change for your own container !
rootmodel=DTINet # change for your own container !

# copy files from machine to docker container
docker cp $folder_path/tmp_data $container:/$rootmodel/
docker cp compute_similarity.m  $container:/$rootmodel/

# once done is, follow the second step as in the readme
