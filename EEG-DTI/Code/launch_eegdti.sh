#!/bin/bash

# launch with nohup 
# nohup bash ./launch_eegdti.sh -b NR -d eeg_dti > launch_eegdti_nr.out &
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


# echo $wdir
folder_path=../Data/$wdir
echo 'directory path: ' $folder_path
 

# 1.Create the container from the DTI2Vec image
eval "DOCKER_ID=$( docker run -d -t $dockerName )";
echo $DOCKER_ID

# 2. Create folders & copy data to execute the model
# first copy from outside everything that we need

echo 'Copying data folder to container...'
docker cp $folder_path $DOCKER_ID:/EEG-DTI 

echo 'Preparing files in docker...'
docker exec $DOCKER_ID mkdir -p /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/oneTooneIndex /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/sevenNets /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/sim_network /EEG-DTI/data_dti


# 3. Execute (outside) once prepared
echo 'Launching eegdti on ' $db
time docker exec $DOCKER_ID  python3 main_modified_eegdti.py $db > ../Results/results_$db.out && docker stop $DOCKER_ID
