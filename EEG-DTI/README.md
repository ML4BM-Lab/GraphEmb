
#  Launching the model in Docker

This section assumes that the preprocessing of the data has been already made.
This means that the Data folders have already been filled.

The following script prepares the docker and launches the model inside the docker.
The results are also writen in Results/ folder as an .out file. 

```
nohup bash ./launch_eegdti.sh -b <db name> -d <image_name> > launch_eegdti_nr.out &
```

***** Test with: 
nohup bash ./launch_eegdti.sh -b NR -d ftest_eegdti > launch_eegdti_nr.out &



All needed data should already be in data_dti. 
The python script alone can be run as:
```
python3 main_modified_eegdti.py <database name folder>
```



# Launching the model with splits in Docker

For this part, again splits shoudld be previously created and in Data/ folder.
These can be created using:

!!!! MODIFY FOR RMSD !!!!
```
python3 generate_onetooneindex_splits.py --dbPath <dbname> --split_type <Sp...> -subsampling -rmsd
```

For launching the model, execute the following shell script. 

```
bash launch_eegdti_splits.sh -b NR -d ftest_eegdti -s Sp (-r for rsmd)
```
test with:
nohup bash launch_eegdti_splits.sh -b NR -d ftest_eegdti -s Sp > log_splits_test.out &

# Note

The docker given includes the modified files in order to launch the model ! 

We needed to modify the following in the default run:
    - main .py
        because: 
            - batch size & number of epochs
    - minibatch.py
        because: folder changed for input new data

We need to modify the following files in the splits run:
    - main .py
    - Minibatch
    - Optimizer (0s)



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

2. Run compute_similarity.m in matlab docker.

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
 -> including the compute_similarity.m script 
 -> including now the transpose matrix of mat_drug_protein.txt (mat_protein_drug.txt)

Folders/files that actually use:
- sevenNets: mat_* from Luo (including transpose)
- sim_network == network from compute similarity in DTINet (compute_similarity.m)
- oneTooneIndex: index of positive and negative edges for train/test


