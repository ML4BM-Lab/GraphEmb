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

# create wdir variable here
# wdir = 

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

#echo $wdir
folder_path=../Data/$wdir
echo 'directory path: ' $folder_path
 

# 1.Create the container from the DTI2Vec image
eval "DOCKER_ID=$( docker run -d -t $dockerName )";
echo $DOCKER_ID

# 2. Create folders & copy data to execute the model
# first copy from outside everything that we need
echo 'Copying files to container...'
#docker cp main_modified_eegdti.py  eeg_dti_nr_last_test:/EEG-DTI
#docker cp minibatch.py  eeg_dti_nr_last_test:/EEG-DTI/decagon/deep
docker cp  main_modified_eegdti.py $DOCKER_ID:/EEG-DTI
docker cp  minibatch.py $DOCKER_ID:/EEG-DTI/decagon/deep
#docker cp ../Data/Yamanashi_et_al_GoldStandard/NR eeg_dti_nr_last_test:/EEG-DTI
# data
docker cp $folder_path $DOCKER_ID:/EEG-DTI 
echo 'Preparing files in docker...'
docker exec $DOCKER_ID mkdir -p /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/oneTooneIndex /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/sevenNets /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/sim_network /EEG-DTI/data_dti


# 3. Execute (outside) once prepared
echo 'Executing eegdti on ' $db
time docker exec $DOCKER_ID  python3 main_modified_eegdti.py $db > $folder_path/results_$db.out && docker stop $DOCKER_ID
