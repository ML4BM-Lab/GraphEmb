#!/bin/bash


while getopts "b:p:d:" opt; do
  case $opt in
    b) db=$OPTARG   ;; # database name
    p) path=$OPTARG   ;; # database name
    d) dockerName=$OPTARG      ;;
    *) echo 'error' >&2
       exit 1
  esac
done

# decide later what I prefer, all path or just DB name
echo "dockerName: $dockerName";
echo "path: $path";
echo "db: $db";


# docker run --name eeg_dti_nr_last_test -it eeg_dti bash
# but with DTI2Vect stuff
# # copy files from machine to docker container

#docker cp main_modified_eegdti.py  eeg_dti_nr_last_test:/EEG-DTI

#docker cp ../Data/Yamanashi_et_al_GoldStandard/NR eeg_dti_nr_last_test:/EEG-DTI
# inside docker create a folder
# mkdir data_dti
# and cp to data_dti the 3 needed folders
# once prepared
# execute outside docker 

time docker exec -it eeg_dti_nr_last_test  python3 main_modified_eegdti.py > results_nr.out
