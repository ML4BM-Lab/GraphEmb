Molecular Docking for subsampling. 
======
# General notes

Joining all proteins by UniprotID, we count **6168** proteins.

From these, we were able to retrieve from the PDB Database and clean **2556** proteins (X-Ray, <2 A ...).
From AlphaFold we take **1777** proteins (restricting to a average of 70 ppXXXX)

This makes that we have structure for a 70.25% of the total number of proteins.

## Codes


## Downloading Data
### PDB Database
Download structures from PDB using the download script.
Before, run get_list of proteins to know which proteins should be download
[batch download](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script)

and download with:

```
nohup bash batch_download.sh -f list_download_pdbs.txt -o PDB_FILES -p > get_pdbs.out &
```



### AlphaFold
All data from AlphaFold downloaded from:
--->  [download page](https://alphafold.ebi.ac.uk/download)
and selected the human prot for later filter.

This wll download everythin, delete *.cif.gz



## Preprocessing data
Select those that we need as in script XXXX

Using Bio3D R package to calculate RMSDs


Docking part? using autodock vina


First we selected and cleaned al available pdbs from PDB database.
Later, we searched the remaining ones in AlphaFold.

Once this is done (all files decompresed in format .pdb)
Presented in folders Clean_from_AFold, Clean_from_PDB


## Calculating RMSD
Using R package bio3d, from it the following funcions (information taken from bio3d
grom the grantlab)

- pdbsplit: 
This function will produce single chain PDB files from multi-chain input files. By default all separate filenames are returned. To return only a subset of select chains the optional input ‘ids’ can be provided to filter the output (e.g. to fetch only chain C, of a PDB object with additional chains A+B ignored). See examples section for further details.

- pdbaln:  
Create multiple sequences alignments from a list of PDB files returning aligned sequence and structure records. 
    * This function uses muscle. relative path can be specified.
    * ncore permits to parallel with 'parallel' installed. 

- rmsd:  Calculate the RMSD between coordinate sets.
    * a numeric vector containing the reference coordinate set for comparison with the coordinates in b. Alternatively, if b=NULL then a can be a matrix or list object containing multiple coordinates for pairwise comparison.
    * fit (logical), if TRUE coordinate superposition is performed on the input structures.
    * ncore: number of CPU cores used to do the calculation. ncore>1 requires package ‘parallel’ installed.


### Download musclev3
[download muscle v3](https://drive5.com/muscle/downloads_v3.htm)