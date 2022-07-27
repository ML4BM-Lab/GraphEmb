
# Instructions for launching the model with docker

This section assumes that the preprocessing of the data has been already made.
This means that the folders in ... are already present

If not, read below.

all needed data should already be in data_dti

```
python3 main_modified_eegdti.py <database name folder>
```

or with docker with the script
```
nohup bash ./launch_eegdti.sh -b NR -d eeg_dti > launch_eegdti_nr.out &
nohup bash ./launch_eegdti.sh -b <database name> -d eeg_dti > launch_eegdti_nr.out &
launch_eegdti.sh <DB name>
```

how all executed:
```
nohup bash ./launch_eegdti.sh -b E -d eeg_dti > log_launch/launch_eegdti_E.out &
nohup bash ./launch_eegdti.sh -b IC -d eeg_dti > log_launch/launch_eegdti_IC.out &
nohup bash ./launch_eegdti.sh -b GPCR -d eeg_dti > log_launch/launch_eegdti_GPCR.out &
nohup bash ./launch_eegdti.sh -b NR -d eeg_dti > log_launch/launch_eegdti_NR.out &

nohup bash ./launch_eegdti.sh -b BindingDB -d eeg_dti > log_launch/launch_eegdti_BindingDB.out &
nohup bash ./launch_eegdti.sh -b Davis_et_al -d eeg_dti > log_launch/launch_eegdti_Davis_et_al.out &
nohup bash ./launch_eegdti.sh -b BIOSNAP -d eeg_dti > log_launch/launch_eegdti_BIOSNAP.out &
nohup bash ./launch_eegdti.sh -b DrugBank -d eeg_dti > log_launch/launch_eegdti_DrugBank.out &
```


Missing:
DrugBank -> launched at  9:44




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
run_preprocess.py <db_name>
```

This model needs 3 folders with specific data:
    - sevenNets
    - sim_network
    - oneTooneIndex 


sim_network corresponds to the output of compute_similarity.m script from DTINet.
In order to retrieve these 3 folders, you should follow some steps, as matlab fails in execution form outside docker 
in some machines. 

1. With this first script we get the folder sevenNetworks, oneTooneIndex, and send data to the dockerfile for executing matlab
```
bash processing_matlab_A.sh <database name>
```

2. Run compute_similarity.m in matlab docker

```
matlab -nodisplay -nosplash -nodesktop -r "run('compute_similarity.m');exit;"
```

This script has been slightly modified for usage issues, only in paths. 
Once run, you can leave the docker and from code folder again you can execute the second
part of the bash file. 


3. Now we have already in the folder the inputs of the model:sevenNets & sim_network with running the following:

```
bash processing_matlab_B.sh <database name>
```


## Data Information
EEG-DTI uses the Luo Dataset (from DTINet).
 -> but including the compute_similarity.m script 
 -> including now the transpose matrix of mat_drug_protein.txt (mat_protein_drug.txt)

Folders/files that actually use:
- sevenNets: mat_* from Luo (including transpose)
- sim_network == network from compute similarity in DTINet (compute_similarity.m)
- oneTooneIndex: index of positive and negative edges for train/test


