#!/bin/bash

while getopts "b:d:" opt; do
  case $opt in
    b) db=$OPTARG   ;; # database name
    *) echo 'error' >&2
       exit 1
  esac
done

# decide later what I prefer, all path or just DB name
if [ $# -eq 0 ]; then
    echo "No DB argument supplied"; exit
fi

echo "db: $db";
# 1.Create the container with

eval "DOCKER_ID=$(docker run -dt dtinet_matlab )";
echo "DOCKER_ID: " $DOCKER_ID

echo "DTINet $db : $DOCKER_ID  " > "log_docker$db.log"

#docker cp  ../Data/$dbname $dockername:/DTINet
#docker cp  Launch_DTINet_matlab.sh $dockername:/DTINet
#docker cp DTINet.m $dockername:/DTINet/src


# creating wdir variable to copy data
yam='Yamanashi_et_al_GoldStandard'
if  [[ "$db" == "NR" ]] || [[ "$db" == "E" ]] || [[ "$db" == "GPCR" ]] || [[ "$db" == "IC" ]]; then
    wdir=$yam/$db
else
    wdir=$db
fi

echo $wdir


# copy files
echo "sending files to docker..."
docker cp  ../Data/$wdir $DOCKER_ID:/DTINet
docker cp  Launch_DTINet_splits_matlab.sh $DOCKER_ID:/DTINet
docker cp DTINet.m $DOCKER_ID:/DTINet/src/DTINet.m

echo "Done! run all in your docker!"

echo "DOCKER_ID: " $DOCKER_ID
