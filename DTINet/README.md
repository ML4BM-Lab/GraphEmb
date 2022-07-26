# Launching DTINet

## Default Settings 

#### Launch DTINet with Docker

Matla gives us problems with x11, 
for that reason the .sh is written to be executed inside the docker. 

The first, you should create a docker:

##### Create container from image (specifying name)
```
docker run -dt --name dtinet_original -it dtinet_matlab bash
```

with dt keeping the docker up, other option is to do:

```
docker restart dtinet_original 
```

#### Copy files
```
docker cp  ../Data/DrugBank dtinet_drugbank:/DTINet
docker cp  Launch_DTINet_matlab.sh dtinet_drugbank:/DTINet
```

#### Enter docker
Enter in the docker interactively and execute the bash script

```
docker exec -it dtinet_original bash
```

##### Run DTINet Matlab Scripts

First time in docker, Matlab will ask for a Email & Pasword

Run specifiying the folder copied as path to copy files
```
bash Launch_DTINet_matlab.sh <DB_folder_name>
```
for example, if we copy the data for DrugBank, we call it as bash Launch_DTINet_matlab.sh DrugBank


This will output a log file: log_DTINet.out
that can be copied to our machine (from outside) as:
```
docker cp dtinet_drugbank:/DTINet/log_DTINet.out <desired_path/log_DTINet_DB.out>
```


### stop docker
If we want to stop our container:
```
docker stop dtinet_drugbank
```




## Evaluation with new splits

### Generating splits
Genereating all splits for a folder can be done with:
```
pyhton3 generate_all_splits_dtine.py <db name>
```

Splits can be generated individually with the script generate generate_splits_dtinet.py.

```
python3 generate_splits_dtinet.py --dbPath <database name> --split_type <Sp/Sd/St>  -subsampling
```

### Launch Model in Docker

Once all splits in folder, you can prepare your docking by runing the following script:

```
bash prepare_launch_splits.sh -b <db name> 
```
this will generate a dockerID were you should execute the folowing in script.
This id is also saved in a logfile. This is because matlab cannot be execute from outside the docker for us.

Then you can enter in the docker with

```
docker exec -it <generated DockerID> bash
```

once in the docker you only have to run the following shell script
specifying the model you want.
Remember that first time will ask for your matlab credentials!

```
bash Launch_DTINet_splits_matlab.sh -b <dbname> -s <Sp/Sd/St> -u <true/false>
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