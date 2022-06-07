

Molecular Docking 
for subsampling approaches. 



## Downloading Data
### PDB Database
Download structures from PDB using the download script.
Before, run get_list of proteins to know which proteins should be download
[batch download](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script)

and download them as

```
nohup bash batch_download.sh -f list_download_pdbs.txt -o PDB_FILES -p > get_pdbs.out &
```



### AlphaFold
All data from AlphaFold downloaded from:
--->  [download page](https://alphafold.ebi.ac.uk/download)
and selected the human prot for later filter.



## Preprocessing data
Select those that we need as in script XXXX

Using Bio3D R package to calculate RMSDs

Docking part? using autodock vina


First we selected and cleaned al available pdbs from PDB database.
Later, we searched the remaining ones in AlphaFold.


## Calculating RMSD

### Download musclev3

[download muscle v3](https://drive5.com/muscle/downloads_v3.htm)