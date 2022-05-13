#!/bin/bash

# to be executed in Code folder! 
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


# copy from docker
container=matlab_env_t2 # change for your own container !
rootmodel=DTINet # change for your own container !

docker cp  $container:/$rootmodel/sim_network $folder_path/
rm -r $folder_path/tmp_data # not needed anymore

# remove files for next use
docker exec $container rm -rf /$rootmodel/compute_similarity.m 
docker exec $container rm -rf /$rootmodel/sevenNets
docker exec $container rm -rf /$rootmodel/tmp_data
docker exec $container rm -rf /$rootmodel/sim_network

docker stop $container

echo "Done!"
echo "You should have now all needed input folders in $wdir"