HyperAttentionDTI
===

# Launch HyperAttentionDTI with Docker

Due to the different libraries this method requires, we are providing a docker to run it over different datasets.

## Prepare docker

Build image as (in the same folder as the Dockerfile):

```
# docker build -t <image_name> .
docker build -t hadti_model .
```

Run the container as:

```
# docker run -dt --name <container_name> <image_name>
docker run -dt --name hadti_benchmark hadti_model
```

## Load Data

The first step is to copy the data files from the database to be used to the docker container. As HyperAttentionDTI only requires the SMILE sequence of the Drug
and the aminoacid sequence of the protein, these have been concatenated into a single file, containing these information plus the interaction label.
For each database, the Sp, Sd, St have been also generated and are located under each database folder. We will use Yamanishi_NR database and Sp mode as an example. These commands should be run from the HyperAttentionDTI's folder.

First, create a folder in the docker's working directory. 

```
mkdir DATASETS
mkdir DATASETS/Yamanishi_NR
mkdir DATASETS/Yamanishi_NR/Sp
```

Outside of the docker, copy files and splits:

```
# copy the splits (this example copies the Sp split)
docker cp Data/Yamanashi_nr/Sp/* hadti_benchmark:/HpyerAttentionDTI/DATASETS/Yamanishi_NR/Sp/.
```

## Launch the model

As this model does not incorporate Sp, Sd and St splits by default, we have modified the main .py file and should be replaced in the docker's container.
It is located under the Code folder, named "HyperAttentionDTI_main.py". This will overwrite the main file of the original HyperAttentionDTI repository.

```
docker cp Code/HyperAttentionDTI_main.py hadti_benchmark:/HpyerAttentionDTI/.
```

To run the model, execute the following command from the folder where the main is located:

```
# lets create a logs and results folders
mkdir Results
mkdir logs

# run the modified version of the NeoDTI main code: 
nohup python -u src/HyperAttentionDTI_main.py Yamanishi_nr Sp 0 > logs/NR_default.out &
```

The first argument specifies the name of the folder containing the database to be evaluated,
the second argument is the mode type (Sp, Sd or St) and the thir argument is the seed numer (0 to 4 both included).
This will run HyperAttentionDTI on Yamanishi NR dataset, split Sp and first seed, creating two folders (<dataset>/<seed>)
where it will save several checkpoint models and results including AUC and AUPRC metrics.

# Code files description

For every database, the aminoacid and SMILES sequence of proteins and drugs, respectively, have been generated. Therefore,
the Code folder of this model contains one script for every database generating the aforementioned inputs.