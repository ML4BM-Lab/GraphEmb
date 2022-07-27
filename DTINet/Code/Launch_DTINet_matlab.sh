# !/bin/bash

DIR_DB=/DTINet/$1 # parse db name

# check data folder
DIR_DATA=/DTINet/data
if [ -d "$DIR_DATA" ];
then
    echo "data folder exist, removing previous files in data"
    rm $DIR_DATA/*.txt
else
    echo "creating data folder"
    mkdir $DIR_DATA
fi

echo "Copy data from $DIR_DB"
cp $DIR_DB/mat*.txt $DIR_DATA/
cp $DIR_DB/Sim*.txt $DIR_DATA/

# initialize network
echo "removing previous files in network and feature"
DIR_NETWORK=/DTINet/network
[ -d "$DIR_NETWORK" ] && rm $DIR_NETWORK/*.txt
DIR_FEAT=/DTINet/feature
[ -d "$DIR_FEAT" ] && rm $DIR_FEAT/*.txt

echo "Launching DTINet..."
echo "Compute similarity networks"
matlab -nodisplay -nosplash -nodesktop -r "run('src/compute_similarity.m');exit;"

echo "Run DCA"
matlab -nodisplay -nosplash -nodesktop -r "run('src/run_DCA.m');exit;"

echo "Launching model with nohup, output in log_DTINet_$1.out"
nohup  matlab -nodisplay -nosplash -nodesktop -r "tic;run('src/run_DTINet.m');toc;exit;" > "log_DTINet_$1.out" 2>/dev/null &
