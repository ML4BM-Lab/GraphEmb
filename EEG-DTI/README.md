

# Instructions for preprocesing data

This model needs the same preprocessing as for DTINet, indeed the XXX Luo Dataset. 

Then we need to apply to the interacction/association/XXXXX matrix, the scripts compute_similarity.m
as in DTINet.

In order to do that, you should follow some steps, as matlab fails in execution form outside docker 
in some machines. 

1. With this first script we get the folder sevenNetworks, and send data to the dockerfile for executing matlab
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

---------------------------------------------------> oneTooneIndex !!!!!!!!!!!!?????


## Data
EEG-DTI uses the Luo Dataset (from DTINet).
 -> but including the compute_similarity.m script 
 -> including the transpose matrix of mat_drug_protein.txt (mat_protein_drug.txt)


=> Any change??? check! 


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