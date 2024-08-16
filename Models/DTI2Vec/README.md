DTI2vec
====

## Prepare docker

Build the image provided in the DockerFile folder:

```
docker build -t <image_name> ./Models/DTI2Vec/Dockerfile/Dockerfile
```
For example: ```docker build -t dti2vec ./Models/DTI2Vec/Dockerfile/Dockerfile.```

Run the container as:

```
docker run -dt --name <container_name> <image_name>
```

For example ```docker run -dt --name dti2vec dti2vec_Sp```.

If you want to run your docker in a given folder from your computer, remember to add the flag ``-v`` indicating the paths to connect.


# Inputs needed for DTI2vec

Raw data needed:
	* Target - Drug interaction list
	* Target - Drug adjacency matrix (called admat)

Processed data:
To create the adjacency matrix, we need to generate the following data based on the Target and Drug properties:

	- Similarity matrix from targets (Smith-Waterman normalized score)
	- Similarity matrix from drugs (SIMCOMP score or similitude matrix between drugs)
	- Drug Target Interaction matrix
	- Drug and Protein embeddings from Node2Vec
 
 The code for creating these matrix for each data base is provided in the Code folder. 
For generating the embeddings, the node2vec software needs to be installed and provided to the corresponding scripts.
 We recommend the node2vec implementation distributed in the SNAP suite.

## Download the data ready to run the models

The first step is downloading the data files from Zenodo (*DTI2Vec_data.tar.gz*).
For each database, both the Sp, Sd, St have been generated and the train/test/val.csv splits.
Those files are located following the structure, using for example Davis and Sp, *Data/Davis/Sp/seed1/train.csv*.

## Run the model inside the docker

```
 docker exec dti2vec python3 -u ./DTi2vec_main_SD_ST_SP.py --data Custom --classifier xgbc --func $func --mode $mode
```
Where the different parameters are:

- ``$func`` is the fusion function and can take two values: Concat or Hadmard
- ``$mode`` corresponds to the splitting mode desired and can take  Sp, St or Sd

## Run the model automatically to reproduce the manuscript results

Finally, we provide a script (```Models/DTI2Vec/Code/Launch_DTI2Vec.sh```) to automate the grid search and comparison between different parameters and databases to be executed on the dockers.

To run this script, the two different matrices need to be generated following this folder order:
```
.
├── Code
├── Data
│   ├── BindingDB
│   │   ├── BindingDB_admat.tsv
│   │   ├── BindingDB_dti.tsv
│   ├── Yamanashi_et_al_GoldStandard
│   │   ├── NR
│   │   │   ├── NR_admat.tsv
│   │   │   ├── NR_dti.tsv
```

Then this can be executed as:
```
nohup bash ./Data/Code/Launch_DTI2Vec.sh -p ./Data/Davis_et_al -d dti2vec:1.0 -m Sp &
nohup bash ./Data/Code/Launch_DTI2Vec.sh -p ./Data/Yamanashi_et_al_GoldStandard/NR -d dti2vec:1.0 -m Sp &
```
