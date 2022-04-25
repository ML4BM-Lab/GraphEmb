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
-- DTI tsv ! 
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

##### Create container from image (specifying name)

```
$ docker run --name dtinet_original -it dtinet_matlab bash
```
.... WORKING HERE ....

##### Keep
```
$ docker restart dtinet_original 
```

#### Copy files
```
docker cp  /home/uveleiro/data/jfuente/DTI/Input4Models/DTINet/Data/DrugBank dtinet_drugbank:/DTINet
```

#### entrar
docker exec -it dtinet_original bash

##### remove files in folders if necesary (data, feature, network)
avoid any problem. safety check,

#### copy files from folder to data
mat*.txt, Sim*.txt

#### Run DTINet Matlab Scripts
First time it will ask for a Email & Pasword
```
matlab -nodisplay -nosplash -nodesktop -r "run('src/compute_similarity.m');exit;"

matlab -nodisplay -nosplash -nodesktop -r "run('src/run_DCA.m');exit;"

matlab -nodisplay -nosplash -nodesktop -r "run('src/run_DTINet.m');exit;"
```
last one mayabe is better to execute as nohup


#### Outside docker & nohup
```
nohup docker exec nombre_container python3 /..../..../.py > path2file.out &
#el file.out se guarda en margaret
```

#### inside docker
nohup  matlab -nodisplay -nosplash -nodesktop -r "run('src/run_DTINet.m');exit;" > log_DTINet_DrugBank.out &
