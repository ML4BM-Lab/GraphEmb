
# Instructions for launching the model with docker

This section assumes that the preprocessing of the data has been already made.
If not, read below.

execute 
```
.sh
```


# Instructions for running the preprocesing of the data

This model needs the same preprocessing as for DTINet, indeed the  Luo Dataset. 

If you want to replicate this step, you can run one of the following scripts
```
***
```

This model needs 3 folders with specific data:
    - sevenNets
    - sim_network
    - oneTooneIndex 


sim_network correspondons to the output of compute_similarity.m script from DTINet.
In order to retrieve these 3 folders, you should follow some steps, as matlab fails in execution form outside docker 
in some machines. 

1. With this first script we get the folder sevenNetworks, oneTooneIndex, and send data to the dockerfile for executing matlab
```
sh processing_matlab_A.sh
```

2. run compute_similarity.m in matlab docker
```
matlab -nodisplay -nosplash -nodesktop -r "run('compute_similarity.m');exit;"
```
This script has been slightly modified for usage issues, only in paths 

3. Now we have already in the folder the inputs of the model:sevenNets & sim_network with running the following:
```
sh processing_matlab_A.sh
```




## Data Information
EEG-DTI uses the Luo Dataset (from DTINet).
 -> but including the compute_similarity.m script 
 -> including now the transpose matrix of mat_drug_protein.txt (mat_protein_drug.txt)


Folders/files that actually use:
- sevenNets: mat_* from Luo (including transpose)
- sim_network == network from compute similarity in DTINet (compute_similarity.m)
- oneTooneIndex:
    * train_index_(0, 1)'+str(seed)+'.txt' for each seed


## Folder structure for input

They have 4 folders in data, from which they only use SevenNets & sim_network, oneToon

## Notes on original dataset & Github

- Yamanishi: error ????  --> in Issues Github! --> ACTUALLY we dont need this one
- Luo: works!

Run without GPU: stopped at seed 10
Run with GPI

## Execute docker

image name: eeg_dti

##### run

ollo ó piollo: run with gpu
pero no nos da la memoria a nosotros, asi que a funcionar sin gpu
```
docker run  --gpus all --name eeg_dti_test_gpu  -it eeg_dti  bash
```

docker run  --name eeg_dti_test_index  -it eeg_dti  bash


inside, can be checked as nvidia-smi 


docker run  --name eeg_dti_test_cleanfolder  -it eeg_dti  bash

--
this is 1.15.5 tensorflow
docker run --name eeg_dti_test -it eeg_dti:1.0 bash

or for tensorflow 1.15.0

docker run --name eeg_dti_test_latest -it eeg_dti:latest bash

##### Keep up
```
$ docker restart eeg_dti_test_gpu 
```

#### entrar
```
docker exec -it dtinet_original bash
```



#### execute outside docker
```
nohup docker exec eeg_dti_test python3 /EEG-DTI/main_luo_all_networks.py > log_EEG-DTI_test_noGPU.out &
```


### Copy Dataset Folder 

--> example from margaret to docker

```
docker cp  /home/uveleiro/data/jfuente/DTI/Input4Models/DTINet/Data/Yamanashi_et_al_GoldStandard/NR eeg_dti_test:/EEG-DTI
```


--> In EEG-DTI we need the similarity networks from DTINet, -------------------------------> **** 


my new test in
docker run  --name eeg_dti_test_one  -it eeg_dti  bash

docker cp  /mnt/md0/data/jfuente/DTI/Input4Models/EEG-DTI/Data/Davis_et_al eeg_dti_test_one:/EEG-DTI

error of sizes !! check code!! --> run in error, done for luo with special shape 708 drugs 1512 proteins 

docker cp  eeg_dti_test:/EEG-DTI/main_luo_all_networks.py   .


main_modified_eegdti.py 
is a modification of main_luo_all_networks.py
because we cant used directly their model
modifications in lines
L292 - L311 which are the input sizes from matrix ! 


docker cp  main_modified_eegdti.py eeg_dti_test_one:/EEG-DTI
