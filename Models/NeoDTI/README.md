NeoDTI
====

# Launch NeoDTI with Docker

Due to the different libraries this method requires, we are providing a docker to run it over different datasets.

## Prepare docker

Build image as (in the same folder as the Dockerfile):

```
# docker build -t <image_name> .
docker build -t neodti_model .
```

Run the container as:

```
# docker run -dt --name <container_name> <image_name>
docker run -dt --name neodti_benchmark neodti_model
```

## Load Data

The first step is to copy the data files from the Formatted folder of the database to be used to the docker container. For each database, the Sp, Sd, St have been also generated and are located under each database folder. We will use Yamanishi_NR database as an example. These commands should be run from the Neo DTI's folder.


First, create a folder in the docker's working directory. 

```
mkdir Yamanishi_NR
```

Outside of the docker, copy files and splits:

```
# copy formatted folder
docker cp Data/Yamanashi_et_al_GoldStandard/NR/Formatted/* neodti_benchmark:/NeoDTI/data/Yamanishi_NR/.

# copy also the splits
docker cp Data/Yamanashi_et_al_GoldStandard/NR/nr_splits* neodti_benchmark:/NeoDTI/data/Yamanishi_NR/.
```

## Launch model

As the NeoDTI model does not include Sd, Sp and St splits by default, we provide a new file in which we included them (default mode as well), that can be located under the Code folder. We will copy it to the code folder of the created docker.

```
docker cp Code/neodti_cv_mod.py neodti_benchmark:/NeoDTI/src/.
```

Once that is done, we can run the model:

```
# enter the docker
docker exec -it neodti_benchmark bash

# change directory to the NeoDTI root folder
cd ./../.

# lets create a logs folders
mkdir logs

# run the modified version of the NeoDTI main code:
nohup python -u src/NeoDTI_cv_mod.py -m Yamanishi_et_al_GoldStandard/NR -s default > logs/NR_default.out &
```

Even though we used the default parameters to run NeoDTI for each database, they can be easily modified:

- **Embedding Dimension (-d or --d)**: 
  - **Default**: 1024 
  - **Description**: Specifies the embedding dimension `d`.

- **Global Norm (-n or --n)**: 
  - **Default**: 1 
  - **Description**: Specifies the global norm to be clipped.

- **Project Matrices Dimension (-k or --k)**: 
  - **Default**: 512 
  - **Description**: Specifies the dimension of project matrices `k`.

- **Test Scenario (-t or --t)**: 
  - **Default**: "o" 
  - **Description**: Specifies the test scenario.

- **Positive Negative Ratio (-r or --r)**: 
  - **Default**: "ten" 
  - **Description**: Specifies the positive-negative ratio.

- **Selected Dataset Matrix (-m or --m)**: 
  - **Default**: "default_data" 
  - **Description**: Specifies the selected dataset matrix. Select the name of the folder within "data" to run NeoDTI on a specific database.

- **Type of Split (-s or --s)**: 
  - **Default**: "default" 
  - **Description**: Specifies the type of split. Set "sp", "sd" or "st" to use each of the Sp, Sd or St splits. Otherwise the default (random) will be used.

# Code files description

This model requires several side-matrices coming from drug-drug interactions, protein-protein interactions, etc. Therefore, many different methods have been run to generate multiple matrices. Also, it has required preprocessing and postprocessing scripts in order to include MOL files of proteins or ensure matrices have a specific format and same number of nodes accross side matrices. Here you can find a description of all the models that have been run and the scripts utilized along the way. 

## Drug Kernels
1. **Drugs_kernels_MorganFingerprint.py**
   - Script for generating drug kernels using Morgan fingerprints.

2. **Drugs_kernels_SIDER.py**
   - Script for generating drug kernels from SIDER data.

## Protein Kernels
1. **Protein_kernels_SW.py**
   - Script for generating protein kernels using Smith-Waterman alignment scores.

## Splits
1. **generate_splits.py**
   - Script for generating Sp, Sd and St splits. They are already provided in every database folder.

## Model
1. **generate_side_matrices.py**
   - Main script to generate all the required side-matrices.
2. **neodti_cv_mod.py**
   - Modified script for cross-validation in NeoDTI models, including Sp, Sd and St splits. Should be copied to docker's container prior to running NeoDTI (see above). 

## Preprocessing & Postprocessing
1. **Postprocessing_format_matrices.py**
   - Script for formatting matrices to a standard structure required for further analysis.
2. **Preprocessing_genMOLfiles.py**
   - Script for preprocessing and generating molecular files necessary for analysis.

## Miscellaneous
1. **Evaluation_augmented_datasets.py**
   - Script for measuring statistics of augmented networks (lost interactions, lost nodes, etc...)


