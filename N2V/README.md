# Node2vec+NN Baseline

# Launch Node2vec+NN with Docker

Due to the different libraries this method requires, we are providing a docker to run it.

## Prepare docker

Build image as (in the same folder as the Dockerfile):

```
# docker build -t <image_name> .
docker build -t N2V_model .
```

Its desired to bind our N2V folder to the docker container, we would use the -v option:

```
# docker run -dt -v <path_to_folder>:/wdir/ --name <container_name> <image_name>
docker run -dt -v <your_path>/N2V/:/wdir/ --name N2V_evaluation N2V_model
```

Once the container is ready, please make sure you have downloaded from Zenodo the data folder
of N2V, and place it within the data folder of N2V. Specifically, the folders "DBs_dict" and
"DBs_embeddings" containing the DTI networks and the output of N2V embedding method for each dataset. There will be many files
as a grid search was run for N2V. Lets now run the provided .py files in the code folder:

## Launch the classification model

First, run "apply_gridsearch.py", where several configurations of N2V will be tried, 
such as the embedding size of the number of steps in the random walk. Also, a few parameters of the
Shallow Neural Network will be included in the grid search, such as the loss type or the type of architecture.
This will generate a large list of pickle files in the Results folder, containing evaluation metrics of the classification 
task for each of the databases, such as AUC or AUPRC. This will take a while so we can save the log in an output file:

```
cd <your_path>/N2V/.
mkdir Results/n2v_nn_results/
python3 -u apply_grid.search.py > Results/grid_search.out &
```

Alternatively, **these files are already provided in the n2v_nn_results folder of the associated Zenodo repository.**
Then, run "group_results.py", that will generate a .tsv file containing the evaluation metrics and the run information 
(N2V parameters) from the generated list of pickle files.

```
python3 -u group_results.py
```

Finally, run "performing_sh_all.py" where from the previous .tsv file, it will pick the best models
from the validation of a network, and it will compute the perfomance of the classification task on the rest of
the networks, generating the Sh matrix.

```
python3 -u performing_sh_all.py
```