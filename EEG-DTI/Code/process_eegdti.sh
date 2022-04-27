#!/bin/bash          
wdir=Davis_et_al # define outside for each database 

# copy necesary matrix in datafolder
folder_path=../Data/$wdir

######### sevenNets  ########
[ ! -d $folder_path/sevenNets ] && mkdir $folder_path/sevenNets
cp $folder_path/mat_*.txt $folder_path/sevenNets

# data for similarity == Luo compute sim, retrieve only \network
[ ! -d $folder_path/tmp_data ] && mkdir $folder_path/tmp_data

cp $folder_path/mat_*.txt $folder_path/tmp_data
cp $folder_path/Sim*.txt $folder_path/tmp_data


######### sim_network ######
# /network cp as sim_network

# docker 
#### run a matlab is DTINEt
#docker run    -it eeg_dti  # --name matlab_in_eeg_dti
#docker exec matlab_in_eeg_dti bash < 'ss' 

# docker run --name matlab_env_test -it dtinet_matlab bash
# docker restart matlab_env_test
path_tmp_folder=$folder_path/tmp_data
name=matlab_env_test
script=compute_similarity.m
# now directory does now work .. copy matlabenv test 
docker cp  $script $name:/EEG-DTI
docker cp  aux.sh $name:/EEG-DTI

docker cp  $path_tmp_folder $name:/EEG-DTI

docker exec -ti $name sh -c "ls" #check that are in folder
echo '--'
docker exec -ti $name sh -c "ls ./tmp_data" #check that are in folder
echo '--'
echo '-> All stuff needed in folder <--'
docker exec -ti $name sh -c "bash aux.sh" #check that are in folder


pwd

# salir del docker?



####### oneTooneIndex ????????????