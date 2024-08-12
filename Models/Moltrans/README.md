Moltrans
====




## Prepare docker

Build image as (in the same folder as the Dockerfile):

```
docker build -t <image_name> .
```
For example: ```docker build -t moltrans .```.

Run the container as:

```
docker run -dt --name <container_name> <image_name>
```

For example ```docker run -dt --name moltrans moltrans_Sp```.

If you want to run your docker in a given folder from your computer, remember to add the flag ``-v``indicating the paths to connect.

## Dowload the data

The first step is to dowload the data files from Zenodo (*Moltrans_data.tar.gz*).
For each database, both the Sp, Sd, St have been generated and the train/test/val.csv splits.
Those files are located following the structure, using as example Davis and Sp, *Data/Davis/Sp/seed1/train.csv*.
Five seeds are available, the same used in the article.
If you are going to 


## Prepare the docker and run the model
If you are running the model inside docker without the ```-v```option.

Outside of the docker, copy data files:

For example, let us say that we want copy the whole data folder, then run the following line:

```
# copy the splits (this example copies the Sp split)
docker cp -r Data/ moltrans:/MolTrans/
```

If we wanted to copy only a given dataset, e.g., Yamanishi NR, then we would need to:

```
# create the folders inside docker
mkdir Data/
mkdir Data/Yamanishi_nr

# copy from outside docker 
docker cp Data/Yamanashi_nr/* moltrans:/MolTrans/Data/Yamanishi_nr/.
```

For running the model with our data, we also need to copy the modified train.py code. 
```
docker cp Code/train.py moltrans:/MolTrans/
``

This version modifies the *get_task* function to enable loading the data from our datafolders. 

To run it, run the following code (following the NR example):


```
python3 train.py --task  'yamanishi_nr_Sp_seed0'
```

In general, for selecting the task or task use the follosing format: <dataset_name> + _ + <split> + _ + seed + <number_seed>, where you can modify the <> fields with your preferences.

Remember that with '-b', '--batch-size' you can specifiy the batchsize (this may be relevant depending on the GPU available).