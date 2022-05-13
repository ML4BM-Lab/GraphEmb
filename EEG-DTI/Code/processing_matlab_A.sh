#!/bin/bash          
# to be executed in Code folder! 

# this XXX for XXX
# change for your XXX
# copy into your machine
# docker cp /home/uveleiro/data/jfuente/DTI/Input4Models/EEG-DTI/Code/process_eegdti_matlab.sh matlab_env_t2:/DTINet


#wdir=$1 # define for each database -> parse this information
db=$1
yam='Yamanashi_et_al_GoldStandard'

if  [[ "$db" == "NR" ]]; then
    wdir=$yam/$db
elif [[ "$db" == "E" ]]; then
    wdir=$yam/$db
elif [[ "$db" == "GPCR" ]]; then
    wdir=$yam/$db
elif [[ "$db" == "IC" ]]; then
    wdir=$yam/$db
else
    wdir=$db
fi

echo $wdir
folder_path=../Data/$wdir

######### sevenNets  #########
[ ! -d $folder_path/sevenNets ] && mkdir $folder_path/sevenNets
cp $folder_path/mat_*.txt $folder_path/sevenNets

### add here oneTooneIndex
python3 generate_onetooneindex.py $db

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

# remove previous files just in case
docker exec $container rm -rf /$rootmodel/compute_similarity.m 
docker exec $container rm -rf /$rootmodel/sevenNets
docker exec $container rm -rf /$rootmodel/tmp_data
docker exec $container rm -rf /$rootmodel/sim_network

# copy files from machine to docker container
docker cp $folder_path/tmp_data $container:/$rootmodel/
docker cp compute_similarity.m  $container:/$rootmodel/

docker restart $container
# once done is, follow the second step as in the readme
echo "Done..."
echo "now you should execute compute similarity script inside matlab docker"
