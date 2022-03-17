# Script List

- process_SIDER
- process_HPRD
- process_Drugbank
- process_CTD
- get_coord_new.py (run processing?)
- generate_all_matrix.py (similarity too? include exc files)

# Input data needed for DTINet

Processed data needed: (download + preproc )

* *** Similarity_Matrix_Drugs.txt (Fingerprint + Tanimoto) ***  ----------------> get_drug_sim.py ---------->
* drug.txt -------------------------------------------------------------> process_DTI_DrugBank.py
* drug_dict_map.txt	----------------------------------------------------> process_DTI_DrugBank.py
* mat_protein_drug.txt -------------------------------------------------> process_DTI_DrugBank.py
* mat_drug_drug.txt ----------------------------------------------------> process_DTI_DrugBank.py
* mat_drug_protein.txt  (transpose) ------------------------------------> --------------------------> preprocess_data_DTINet_DrugBank.py
* mat_drug_protein_remove_homo.txt (remove homologous above) -----------> --------------------------> preprocess_data_DTINet_DrugBank.py

* protein.txt  ---------------------------------------------------------> process_HPRD_DTINet.py
* protein_dict_map.txt -------------------------------------------------> process_HPRD_DTINet.py
* mat_protein_protein.txt ----------------------------------------------> process_HPRD_DTINet.py
* *** Similarity_Matrix_Proteins.txt (Smith-Waterman score) -> script get_SW_score find) *** 

* disease.txt ---------------------------------------------------------> process_CTD_DTINet.py
* mat_protein_disease.txt  --------------------------------------------> process_CTD_DTINet.py (only coordinates)
* mat_drug_disease.txt	

* se.txt  -------------------------------------------------------------> process_SIDER_DTINet.py
* mat_drug_se.txt  ----------------------------------------------------> process_SIDER_DTINet.py (only coordinates)


Original DBs:
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

