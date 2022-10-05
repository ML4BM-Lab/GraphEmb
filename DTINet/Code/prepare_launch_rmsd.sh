#!/bin/bash

while getopts "b:t:p:" opt; do
  case $opt in
    b) db=$OPTARG   ;; # database name
    t) th=$OPTARG   ;; # database name
    p) pn=$OPTARG   ;; # database name
    *) echo 'error' >&2
       exit 1
  esac
done

# decide later what I prefer, all path or just DB name
if [ $# -eq 0 ]; then
    echo "No DB argument supplied"; exit
fi


echo "db: $db";
echo "threshold: $th"
echo "pos2neg ratio: $pn"


# 1.Create the container 
eval "DOCKER_ID=$(docker run -dt dtinet_matlab )";
echo "DOCKER_ID: " $DOCKER_ID
#echo "DTINet $db : $DOCKER_ID  " > "log_docker_rmsd_$db.out"


# 0. generate splits
# testing splits

# generate_subsampling_dtinet.pt
python3 -u generate_subsampling_dtinet.py --dbPath $db -t $th -p $pn  > log_biosnap_"$th"_"$pn"_"$DOCKER_ID".out



# 2. Copy files to docker:

# Creating wdir variable to copy data
yam='Yamanashi_et_al_GoldStandard'
if  [[ "$db" == "NR" ]] || [[ "$db" == "E" ]] || [[ "$db" == "GPCR" ]] || [[ "$db" == "IC" ]]; then
    wdir=$yam/$db
else
    wdir=$db
fi

echo "wdir: " $wdir

# Copy files
echo "Copying files to docker..."

# Copying the full folder with all splits ! 
docker cp  ../Data/$wdir $DOCKER_ID:/DTINet
# Launcher and code
docker cp  Launch_DTINet_RMSD_matlab.sh $DOCKER_ID:/DTINet
docker cp DTINet_rmsd.m $DOCKER_ID:/DTINet/src/DTINet.m
docker cp run_DTINet_rmsd.m $DOCKER_ID:/DTINet/src/run_DTINet.m

echo "Done!"
echo "Entering in your docker container..."

echo "DOCKER_ID: " $DOCKER_ID

docker exec -it $DOCKER_ID bash
