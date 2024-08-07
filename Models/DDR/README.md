DDR
======

# Launching DDR
## Default Settings 
### Launch DDR with Docker

Due to the different libraries this method requires, we are providing a docker to run it over different datasets.

#### Prepare docker

Build image as (in the same folder as the Dockerfile):

```
# docker build -t <image_name> .
docker build -t ddr_model .
```

Run the container as:

```
# docker run -dt --name <container_name> <image_name>
docker run -dt --name ddr_benchmark ddr_model
```
#### Load Data

The first step is to copy the data files from the PostSNF folder of the database to be used
to the docker container. We will use Yamanishi_NR database as an example. These commands should be
run from the DDR's folder. First create a folder in the docker's working directory.

```
docker exec -it ddr_benchmark bash
mkdir DATASETS
mkdir DATASETS/Yamanishi_NR
```

Outside of the docker, copy the files to be used:

```
docker cp Data/Yamanashi_et_al_GoldStandard/NR/Formatted/PostSNF ddr_benchmark:/ddr/DATASETS/Yamanishi_NR
```

#### Launch model

We modified the original code in DDR to print results, and set parameters as input instead of 
hardcoding them as original described in the code. Therefore, the main file "nr_demo.py" should be replaced
before running the model. It ise available at the code folder, so we should just copy it
into the docker's folder.

```
docker cp Code/nr_demo.py ddr_benchmark:/ddr/Similarity_selection_code/DDR_demo/.
```

Once that is done, we can run the model:

```
docker exec -it ddr_benchmark bash
cd DATASETS/Yamanishi_NR/PostSNF/
nohup python -u ../../../Similarity_selection_code/DDR_demo/nr_demo.py nr_admat_dgc_mat_2_line.txt selected_drugs_post.txt selected_prots_post.txt > results_nr.out &
```

We can move the results file out of the docker into the results folder. Outside of the docker,
we would run the following command:

```
docker cp ddr_benchmark:/ddr/DATASETS/Yamanishi_NR/PostSNF/results_nr.out 
```

## Evaluation with Sp,Sd and St splits

This model already incorporates this splits into its pipeline, so it has not required further modifications.

# Code files description

This model requires several side-matrices coming from drug-drug interactions, protein-protein interactions, etc. Therefore, many different methods have been run to generate multiple matrices. Also, it has required both preprocessing and postprocessing scripts, as it uses the SNF method to weight the input matrices before running the model. This requires the matrices should have a specific format and common nodes, promoting many of the side-matrices to be discarded prior to run the model. Here you can find a description of all the models that have been run and the scripts utilized along the way. 


# Script Descriptions

## Preprocessing
1. **Preprocessing_genMOLfiles.py**
   - Script for preprocessing and generating molecular files necessary for smiles and other side matrices.
   
2. **Preprocessing_genProtSymbols.py**
   - Script for preprocessing and generating protein symbols for the datasets.

## Postprocessing
1. **Postprocessing_check_matrices.py**
   - Script to check the consistency and correctness of the generated matrices.
   
2. **Postprocessing_drugs_AERS.py**
   - Script for postprocessing drug-related data, specifically from the AERS database.
   
3. **Postprocessing_format_matrices.py**
   - Script for formatting matrices to a standard structure required to run DDR.

## Protein Kernels
1. **Protein_kernels_BioGRID.py**
   - Script for generating protein kernels using data from BioGRID.
   
2. **Protein_kernels_bioMART.R**
   - R script for generating protein kernels using data from bioMART.
   
3. **Protein_kernels_helpers.py**
   - Helper functions for processing and generating protein kernels.
   
4. **Protein_kernels_kebabs.R**
   - R script for generating protein kernels using the 'kebabs' package.
   
5. **Protein_kernels_SW.py**
   - Script for generating protein kernels using Smith-Waterman alignment scores.

## Drug Kernels
1. **Drugs_kernels_AERS_FDA.py**
   - Script for generating drug kernels from FDA's AERS data.
   
2. **Drugs_kernels_helpers.py**
   - Helper functions for processing and generating drug kernels.
   
3. **Drugs_kernels_Rchemcpp.R**
   - R script for generating drug kernels using the Rchemcpp package.
   
4. **Drugs_kernels_SIDER.py**
   - Script for generating drug kernels from SIDER data.
   
5. **Drugs_kernels_SIMCOMP.py**
   - Script for generating drug kernels using SIMCOMP data.

## Other Scripts
1. **Evaluation_augmented_datasets.py**
   - Script for evaluating augmented datasets.
   
2. **generate_side_matrices.py**
   - Main script to generate all matrices required for running DDR.
   
3. **nr_demo.py**
   - Main file to run ddr within the docker file. Should be copy and replace inside the docker's folder.


# Entity Kernels Information Source

The total number of kernels used for drugs and proteins is detailed in the following table:

| **Entity** | **Kernel**                        | **Information Source**                       |
|------------|-----------------------------------|----------------------------------------------|
| **Drugs**  | AERS-bit                          | AERS bit Side-effects                        |
|            | AERS-freq                         | AERS freq Side-effects                       |
|            | GIP                               | Gaussian Interaction Profile Network         |
|            | LAMBDA                            | Lambda-k Kernel Chem. Struct.                |
|            | MARG                              | Marginalized Kernel Chem. Struct.            |
|            | MINMAX                            | MinMax Kernel Chem. Struct.                  |
|            | SIMCOMP                           | Graph kernel Chem. Struct.                   |
|            | SIDER                             | Side-effects Similarity Side-effects         |
|            | SPEC                              | Spectrum Kernel Chem. Struct.                |
|            | TAN                               | Tanimoto Kernel Chem. Struct.                |
| **Proteins**| GIP                              | Gaussian Interaction Profile Network         |
|            | GO                                | Gene Ontology Semantic Similarity Func. Annot. |
|            | MIS-k3m1                          | Mismatch kernel (k = 3, m = 1) Sequences     |
|            | MIS-k4m1                          | Mismatch kernel (k = 4, m = 1) Sequences     |
|            | MIS-k3m2                          | Mismatch kernel (k = 3, m = 2) Sequences     |
|            | MIS-k4m2                          | Mismatch kernel (k = 3, m = 2) Sequences     |
|            | PPI                               | Proximity in protein-protein network         |
|            |                                   | Protein-protein Interactions                 |
|            | SPEC-k3                           | Spectrum kernel (k = 3) Sequences            |
|            | SPEC-k4                           | Spectrum kernel (k = 4) Sequences            |
|            | SW                                | Smith-Waterman alignment score Sequences     |
