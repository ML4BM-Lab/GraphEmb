# Preprocessing of Databases

Run (file not available for all databases yet)
```
python3 run_preprocess_Yamanisihi.py 
python3 run_preprocess_DrugBank.py *
python3 run_preprocess_BIOSNAP.py *
python3 run_preprocess_Davis.py *
python3 run_preprocess_BindingDB.py *
```

Script that calls:
1. get_coord.py
    This script obtaind all necesary edge lists from complementary databases
    * Calls:
      1. process_HPRD_DTINet.py
      2. process_SIDER_DTINet.py
      3. process_CTD_DTINet.py
      4. process_DrugBank.py
2. get_all_matrix_DrugBank.py
    This script obtains the DTIs and uses the output from previous script to create the needed matrix. 



## List of scripts in 'DTINet/Code/':

- ... FIRST 
- process_SIDER_DTINet.py
- process_HPRD_DTINet.py
- process_Drugbank_DTINet.py
- process_CTD_DTINet.py
- get_coord.py 
- get_all_matrix_DrugBank.py 
- get_all_matrix_Yamanishi.py
- ...

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

## Execute docker


#### Create container from image (specifying name)
```
$ docker run --name dtinet_original -it dtinet_matlab bash
```
.... WORKING HERE ....