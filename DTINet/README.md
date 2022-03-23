# Preprocessing of Databases

preprocessing_DTINet_DrugBank.py (for example, not available yet)

Script that calls:
    1. get_coord.py
        - meaning:
        - Calls:
            1.1. process_HPRD_DTINet.py
            1.2. process_SIDER_DTINet.py
            1.3. process_CTD_DTINet.py
            1.4. process_DrugBank.py
    2. get_all_matrix_DrugBank.py


## Script List

- process_SIDER
- process_HPRD
- process_Drugbank
- process_CTD
- get_coord.py 
- get_all_matrix.py 


# Input data needed for DTINet

Processed data needed: (download + preproc )


## What do we retrieve from databases? 
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

