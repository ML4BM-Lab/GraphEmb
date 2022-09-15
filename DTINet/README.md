# Launching DTINet

## Default Settings 

### Launch DTINet with Docker

Matlab problems with x11 for our workstation. 
For that reason we use 2 .sh codes to prepare the docker and 
to run it for a given dataset.

#### Prepare docker

This script will copy all files and run iteratively the docker.
```
bash ./prepare_dtinet_matlab.sh -b <database_name>  -d <docker_image>
```
in our case docker image is docker_image = dtinet_matlab.

#### Launch model

The first step already copied a .sh launcher, then we need to run it, 
parsing again the db name. 
With this, we could have copied all above instead of only one folder.

It may be needed to run first matlab to activate the software 
in the docker with account e-mail and pasword.

```
nohup bash ./Launch_DTINet_matlab.sh <database name> &
```

for example, if we copy the data for DrugBank, we call it as bash Launch_DTINet_matlab.sh DrugBank.

### Stop docker
If we want to stop our container:
```
docker stop <docker name/docker ID>
```




## Evaluation with new splits


The procedure is similar as before. 

Still, splits can be generated individually with the script generate generate_splits_dtinet.py
This is valid for rmsd as well! it only needs to be parsed

```
python3 generate_splits_dtinet.py --dbPath <database name> --split_type <Sp/Sd/St>  -subsampling -rmsd
```

### Prepare Model 

This split alreay enters in the docker but also saves the docker id.

```
bash prepare_launch_splits.sh -b <db_name> 
```
### Launch Model in Docker

Remember that first time will ask for your matlab credentials!

It can be also launched with nohup ! 

```
bash Launch_DTINet_splits_matlab.sh -b <db_name> -s <Sp/Sd/St> -r
```
-r is an optional argument, if selected the model is executed with rmsd
if not, the normal subsampling is aplied. 


A log file is saved wih a name that specifies the run settings. 


## Testing RMSD

created new code for splits in
```
generate_splits_dtinet.py
```

Then, create a prepare_launch_rmsd.sh
similar to prepare_launch_rmsd.sh
but only call rmsd 

inside this script it will copy the launcher, that also needs to be modify
to copy only the RMSD folder into the splits folder!

then run:
```
bash prepare_launch_rmsd.sh -b BIOSNAP
```
That copy the file Launch_DTINet_RMSD_matlab.sh

```
bash Launch_DTINet_RMSD_matlab.sh
```


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
    - *** DTI *** - change DB in each folder
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


## List of scripts in 'DTINet/Code/':

[]