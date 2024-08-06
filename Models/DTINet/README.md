DTINet
======

# Launching DTINet

## Default Settings 

### Launch DTINet with Docker

Matlab raises problems with x11 for our workstation. 
For this reason, we use two bash codes to prepare the docker and run it for a given dataset.


#### Prepare docker

This script will copy all files and run iteratively the docker.
```
bash ./prepare_dtinet_matlab.sh -b <database_name>  -d <docker_image>
```

#### Launch model

The first step already copied a .sh launcher.
Then we need to run it, parsing again the db name. 

It may be needed to run first matlab to activate the software 
in the docker with account e-mail and pasword.

```
nohup bash ./Launch_DTINet_matlab.sh <database name> &
```

For example:
```
nohup bash ./Launch_DTINet_matlab.sh DrugBank &
```



## Evaluation with new splits

The procedure is similar to the previous one. 

Still, splits can be generated individually with the script generate generate_splits_dtinet.py

```
python3 generate_splits_dtinet.py --dbPath <database name> --split_type <Sp/Sd/St>  -subsampling 
```

Note that the splits were generated using the code of GraphGUEST before it was available as a pypi package.
For using the splits in a new model, it might be faster, simpler and safer to install GUEST (repo: https://github.com/ubioinformat/GraphGuest)


### Prepare Model 

This script generates the splits and a docker container. 
In this step the splits are copied into the docker, and the docker id will also appear on the screen.


```
bash prepare_launch_splits.sh -b <db_name> 
```

### Launch Model in Docker

Remember that first time will ask for your Matlab credentials.

```
bash Launch_DTINet_splits_matlab.sh -b <db_name> -s <Sp/Sd/St> 
```


A log file is saved wih a name that specifies the run settings. 



# Preprocessing of Databases

This model needs different databases. 
All data can be produced running each model independantly or sequentially with the script run_preprocesss_{DB}.py,
where DB = {NR, E, IC, GPCR, DrugBank, BIOSNAP, Davis, BindingDB}.

```
python3 run_preprocess.py <DB>
```

Script that calls:
1. get_coord.py
    This script obtaind all necesary edge lists from complementary databases
    * Calls:
      1. process_HPRD_DTINet.py
      2. process_SIDER_DTINet.py
      3. process_CTD_DTINet.py
      4. process_DrugBank.py 

2. Get DTI (only for Yamanishi, BindingDB, Davis)
    - DTI_Yamanishi.py
    - DTI_BindingDB.py
    - DTI_Davis.py

3. get_all_matrix.py
    This script obtains the DTIs and uses the output from previous script to create the needed matrix. 


## What do we retrieve from each database? 
* DrugBank database 
    - Drug nodes
    - *** DTI *** - change for DB in each folder
    - drug-drug 

* HPRD database
    - Protein nodes
    - Protein-Protein interactions

* Comparative Toxicogenomics Database
    - Disease nodes
    - drug-disease associations
    - protein-disease associations 

* SIDER database.
    - Side-effect nodes
    - drug-side-effect associations 


