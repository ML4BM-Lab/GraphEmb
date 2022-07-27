#!/bin/bash

# use as:
# bash ./prepare_dtinet_matlab.sh -b <DB name> -d dtinet_matlab 

while getopts "b:d:" opt; do
  case $opt in
    b) db=$OPTARG   ;; # database name
    d) dockerName=$OPTARG      ;;
    *) echo 'error' >&2
       exit 1
  esac
done


# decide later what I prefer, all path or just DB name
echo "dockerName: $dockerName";
echo "db: $db";

# creating wdir variable to copy data
yam='Yamanashi_et_al_GoldStandard'
if  [[ "$db" == "NR" ]] || [[ "$db" == "E" ]] || [[ "$db" == "GPCR" ]] || [[ "$db" == "IC" ]]; then
    wdir=$yam/$db
else
    wdir=$db
fi

folder_path=../Data/$wdir
echo 'directory path: ' $folder_path

# 1.Create the container from the matlab image
eval "DOCKER_ID=$( docker run -d -t $dockerName )";
echo $DOCKER_ID > "ID_DOCKER_${db}.out"


# 2. Create folders & copy data to execute the model
# first copy from outside everything that we need
echo 'Copying files to container...'

# Launcher Code
docker cp  Launch_DTINet_matlab.sh $DOCKER_ID:/DTINet
# Data Folder
docker cp $folder_path $DOCKER_ID:/DTINet



docker exec -it $DOCKER_ID bash