# !/bin/bash
## chang this

DIR_DB=/DTINet/$1
# later parse subsampling type 

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

if [ -d splits ];
then
    echo "splits folder exist, removing previous files in data"
    rm splits/*.txt
else
    echo "creating splits folder"
    mkdir splits
fi

echo "Copy data from $DIR_DB"
cp $DIR_DB/mat*.txt $DIR_DATA/
cp $DIR_DB/Sim*.txt $DIR_DATA/
#cp  $DIR_DB/Index_Sp_nosubsampl/*.txt splits/
#cp  $DIR_DB/Index_Sp_subsampl/*.txt splits/
#cp  $DIR_DB/Index_Sd_subsampl/*.txt splits/
cp  $DIR_DB/Index_St_subsampl/*.txt splits/

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

echo "Launching model with nohup, output in 'log_DTINet.out'"
#nohup matlab -nodisplay -nosplash -nodesktop -r "run('src/run_DTINet.m');exit;" > log_DTINet_subs_Sd.out 2>/dev/null &
