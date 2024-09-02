Notes on the calculation of the pairwise RMSD for subsampling. 
======
# General notes

Joining all proteins by UniprotID, we count **6168** proteins.

From there, we could retrieve from the PDB Database and clean **2556** proteins (X-Ray, <2 A).
From AlphaFold we take **1776** proteins (restricting to an average of 70 pLDDT)

Then, we have structure for 70.25% of the proteins in all datasets.

For running these scripts the packages *pymol* and *multiprocessing* are needed!

There are several codes in this folder that were only used for retrieving statistics and ploting data.

## Codes
#### RMSD Calculation

Prioir to the RMSD calculation, the cleaned structures from PDB and AlphaFold simulations are neeeded.

For the RMSD calculation we make use of a modified version of align_all_to_all.py from PyMol scripts.
In this script we consider the superimposition align method.
Changes in the original script are only for saving the output matrix. 


```
calculate_rmsd_pymol.py
```


For the selected proteins, the following script was used to retrieve annotation. 

```
get_annotations_uniprot.py 
```

Pymol may be installed with conda, however there is also availble the pymol-open-source package.


## Downloading Data

### PDB Database

Download structures from PDB using the download script.
Before, run get_list of proteins to know which proteins should be download using
[batch download](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script)

This can be done with with:

```
nohup bash batch_download.sh -f list_download_pdbs.txt -o PDB_FILES -p > get_pdbs.out &
```


### AlphaFold

All data from AlphaFold downloaded from the [download page](https://alphafold.ebi.ac.uk/download)
and selected the humanan proteins for filtering among all the avaliable proteins.

This will perform a full download, however the  *.cif.gz files can be deleted, as we are not using them in our preprocesing pipeline.


## Preprocessing data

Inspect list of proteins and download them from the PDB Database.

```
get_list_proteins_PDB.py
```

From those, retrieve from PDB and clean them, move available files to the Clean_PDB folder.
First we selected and cleaned al available pdbs from PDB database.

```
prepare_PDB_files.py
```

Later, we searched the remaining ones in AlphaFold (the data needs to be bulk downloaded).
These were copied to our work folder and uncompressed. 
```
get_list_proteins_afold.py
```

Once this is done (all files decompresed in format .pdb) structures were saved in in folders Clean_from_AFold, Clean_from_PDB.
And all preprocessing for targets is done. 
