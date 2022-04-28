#!/bin/bash

# to be executed in Code folder! 
wdir=Davis_et_al # define for each database 
folder_path=../Data/$wdir

# copy from docker
container=matlab_env_t2 # change for your own container !
rootmodel=DTINet # change for your own container !

docker cp  $container:/$rootmodel/sim_network $folder_path/
rm -r $folder_path/tmp_data # not needed anymore

# oneToone?