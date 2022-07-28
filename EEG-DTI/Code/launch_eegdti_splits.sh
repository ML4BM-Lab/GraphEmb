#!/bin/bash

# launch with nohup 
# nohup bash ./launch_eegdti.sh -b NR -d eeg_dti -> launch_eegdti_nr.out &
# subsampling is as default

RMSD=false # initialise
while getopts "b:d:s:r" opt; do
  case $opt in
    b) db=$OPTARG   ;; # database name
    d) dockerName=$OPTARG      ;;
    s) SPLIT=$OPTARG      ;;
    r) RMSD=true      ;;
    *) echo 'error' >&2
       exit 1
  esac
done

# decide later what I prefer, all path or just DB name
echo "dockerName: $dockerName";
echo "Database: $db";
# say split type
if [[ "$SPLIT" == 'Sp' || "$SPLIT" == 'Sd' || "$SPLIT" == 'St' ]]; then
    echo "Split type:" $SPLIT
else
    echo "Split not specified"
    exit
fi
# say if rmsd
if "$RMSD"; then
    echo 'Using RMSD'
else
    echo 'Not using RMSD'
fi


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
 

# 1.Create the container from the selected image
eval "DOCKER_ID=$( docker run -d -t $dockerName )";
echo $DOCKER_ID

# 2. Create folders & copy data to execute the model
# first copy from outside everything that we need

echo 'Copying data folder to container...'
docker cp $folder_path $DOCKER_ID:/EEG-DTI 

echo 'Preparing files in docker...'
docker exec $DOCKER_ID mkdir -p /EEG-DTI/data_dti
#docker exec $DOCKER_ID cp -r /EEG-DTI/$db/oneTooneIndex /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/sevenNets /EEG-DTI/data_dti
docker exec $DOCKER_ID cp -r /EEG-DTI/$db/sim_network /EEG-DTI/data_dti

# now copy in one to index the selected splitype
# Copying splits to split folder 
if "$RMSD"; then
  folder_index=oneTooneIndex_"$SPLIT""_subsampling_rmsd"
  docker exec $DOCKER_ID cp -r /EEG-DTI/$db/$folder_index /EEG-DTI/data_dti/oneTooneIndex
else
  folder_index=oneTooneIndex_"$SPLIT""_subsampling"
  docker exec $DOCKER_ID cp -r /EEG-DTI/$db/$folder_index /EEG-DTI/data_dti/oneTooneIndex
fi
echo "copied in OneToOneIndex: $folder_index"


# 3. Execute (outside) once prepared

# modify the log out
if "$RMSD"; then
    out_file=log_EEG-DTI_"$DB"_"$SPLIT"_"subsampling_rmsd.out"
else
    out_file=log_EEG-DTI_"$DB"_"$SPLIT"_"subsampling.out"
fi

echo 'Launching eegdti on ' $db
time docker exec $DOCKER_ID  python3 main_modified_eegdti_splits.py $db > ../Results/$out_file && docker stop $DOCKER_ID
