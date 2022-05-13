
# Instructions for launching the model with docker

This section assumes that the preprocessing of the data has been already made.
If not, read below.

all needed data should already be in data_dti

```
execute main_modified_eegdti.py <DB name>
```

or with docker with the script
```
nohup bash ./launch_eegdti.sh -b NR -d eeg_dti > launch_eegdti_nr.out &
nohup bash ./launch_eegdti.sh -b <database name> -d eeg_dti > launch_eegdti_nr.out &
launch_eegdti.sh <DB name>
```
# Notes

we needed to modify the following .py
    - main **
        because: 
            - batch size & number of epochs
    - minibatch.py
        because: folder changed for input new data

# Instructions for running the preprocesing of the data

This model needs the same preprocessing as for DTINet, indeed the  Luo Dataset. 
If you want to replicate this step, you can run one of the following scripts
```
run_preprocess_***.py
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
sh processing_matlab_A.sh <database name>
```

2. run compute_similarity.m in matlab docker
```
matlab -nodisplay -nosplash -nodesktop -r "run('compute_similarity.m');exit;"
```
This script has been slightly modified for usage issues, only in paths. 
Once run, you can leave the docker and from code folder again you can execute the second
part of the bash file. 


3. Now we have already in the folder the inputs of the model:sevenNets & sim_network with running the following:
```
sh processing_matlab_B.sh <database name>
```




## Data Information
EEG-DTI uses the Luo Dataset (from DTINet).
 -> but including the compute_similarity.m script 
 -> including now the transpose matrix of mat_drug_protein.txt (mat_protein_drug.txt)

Folders/files that actually use:
- sevenNets: mat_* from Luo (including transpose)
- sim_network == network from compute similarity in DTINet (compute_similarity.m)
- oneTooneIndex: index of positive and negative edges for train/test


## Folder structure for input

They have 4 folders in data, from which they only use SevenNets & sim_network, oneToon


## About docker

image name: eeg_dti

if available gpu 
(pero no nos da la memoria a nosotros, asi que a funcionar sin gpu)
```
docker run  --gpus all --name eeg_dti_test_gpu  -it eeg_dti  bash
```

inside, can be checked if the gpu is activated typing nvidia-smi.

```
docker run  --name eeg_dti_test_index  -it eeg_dti  bash
```

-----------------------------------------------------------------------> from here only notes from tests 





#### execute outside docker
```
nohup docker exec eeg_dti_test python3 /EEG-DTI/main_luo_all_networks.py > log_EEG-DTI_test_noGPU.out &
```


### Copy Dataset Folder 

--> example from margaret to docker

```
docker cp  /home/uveleiro/data/jfuente/DTI/Input4Models/DTINet/Data/Yamanashi_et_al_GoldStandard/NR eeg_dti_test:/EEG-DTI
```



main_modified_eegdti.py 
is a modification of main_luo_all_networks.py
because we cant used directly their model
modifications in lines
L292 - L311 which are the input sizes from matrix ! 


docker cp  main_modified_eegdti.py eeg_dti_test_one:/EEG-DTI
