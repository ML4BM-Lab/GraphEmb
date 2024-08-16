DTiGEMS
=======

## Prepare docker

Build the image provided in the DockerFile folder:

```
docker build -t <image_name> ./Models/DTiGEMS/Dockerfile/Dockerfile
```
For example: ```docker build -t DTiGEMS ./Models/DTiGEMS/Dockerfile/Dockerfile.```

Run the container as:

```
docker run -dt --name <container_name> <image_name>
```

For example ```docker run -dt --name DTiGEMS moltrans_Sp```.

If you want to run your docker in a given folder from your computer, remember to add the flag ``-v`` indicating the paths to connect.


# Inputs needed for DTiGEMS

DTiGEMS needs multiple Target similitude  and Drug similitude matrices:
  * Drug - Drug similitude matrix calculated based on FDA_aers_bit
  * Drug - Drug similitude matrix calculated based on FDA_aers_freq
  * Drug - Drug similitude matrix calculated based on Rchemcpp_lambda
  * Drug - Drug similitude matrix calculated based on Rchemcpp_marginalized
  * Drug - Drug similitude matrix calculated based on Rchemcpp_minmaxTanimoto
  * Drug - Drug similitude matrix calculated based on Rchemcpp_spectrum
  * Drug - Drug similitude matrix calculated based on Rchemcpp_tanimoto
  * Drug - Drug similitude matrix calculated based on SIDER_SideEffect
  * Drug - Drug similitude matrix calculated based on SIMCOMP_scores
  * Target - Target similitude matrix calculated based on GO_PPI
  * Target - Target similitude matrix calculated based on mismatch_kernel_k3m1
  * Target - Target similitude matrix calculated based on mismatch_kernel_k3m2
  * Target - Target similitude matrix calculated based on mismatch_kernel_k4m1
  * Target - Target similitude matrix calculated based on mismatch_kernel_k4m2
  * Target - Target similitude matrix calculated based on SmithWaterman_scores_MinMax
  * Target - Target similitude matrix calculated based on spectrum_kernel_k3
  * Target - Target similitude matrix calculated based on spectrum_kernel_k4
  * An adjacency matrix admat_dgc
  * The Drug - Target interaction list

All these matrices need to be located on the same folder.
 
The code for creating these matrices for each database is provided in the Code folder. 
This model needs the embeddings from Node2Vec that are created inside the docker after the new amplified  adjacency matrix is created.

## Download the data ready to run the models

The first step is downloading the data files from Zenodo (*DTiGEMS_data.tar.gz*).
For each database, both the Sp, Sd, St have been generated and the train/test/val.csv splits.
Those files are located following the structure, using for example Davis and Sp, *Data/Davis/Sp/seed1/train.csv*.

## Run the model
To run this model, we provide two different scripts, one using the default splits and other accepting the proposed Sp, Sd, and St splits.
These scripts are provided as Models/DTiGEMS/Code/Launch_DTiGEMS.sh and  Models/DTiGEMS/Code/Launch_DTiGEMS_splits.sh respectively.

The scripts accept as arguments:
- ``DATA_NAME`` refers to the database and the path where to retrieve the matrices.
- ``SPLIT`` referring to the split to be applied
