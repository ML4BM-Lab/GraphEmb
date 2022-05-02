# Preprocessing of Databases

Run (file not available for all databases yet)
'*' missing scripts

```
python3 run_preprocess_Yamanisihi.py 
python3 run_preprocess_DrugBank.py *
python3 run_preprocess_BIOSNAP.py *
python3 run_preprocess_Davis.py 
python3 run_preprocess_BindingDB.py 
```

Script that calls:
1. get_coord.py
    This script obtaind all necesary edge lists from complementary databases
    * Calls:
      1. process_HPRD_DTINet.py
      2. process_SIDER_DTINet.py
      3. process_CTD_DTINet.py
      4. process_DrugBank.py (for)
-- DTI tsv ! 
2. Get DTI (only for Yamanishi, BindingDB, Davis)
    - DTI_Yamanishi.py
    - DTI_BindingDB.py
    - DTI_Davis.py

3. get_all_matrix_DrugBank.py
    This script obtains the DTIs and uses the output from previous script to create the needed matrix. 


## List of scripts in 'DTINet/Code/':

- helper_functions_dtinet.py
- get_coord.py 
- process_SIDER_DTINet.py
- process_HPRD_DTINet.py
- process_Drugbank_DTINet.py
- process_CTD_DTINet.py

- DTI_BindingDB.py
- DTI_Davis.py
- DTI_Yamanishi.py

- get_all_matrix_DrugBank.py 
- get_all_matrix_Yamanishi.py
- get_all_matrix_Davis.py
- get_all_matrix_BindingDB.py
- get_all_matrix_BIOSNAP.py

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


## Launch DTINet with Docker

Matlab problems with x11, for that reason we give the .sh to be executed inside the docker. 

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