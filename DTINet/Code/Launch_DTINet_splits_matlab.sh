# !/bin/bash


while getopts "b:s:u:" opt; do
  case $opt in
    b) DB=$OPTARG   ;; # database name
    s) SPLIT=$OPTARG      ;;
    u) SUBSAMPL=$OPTARG      ;;
    *) echo 'error' >&2
       exit 1
  esac
done

# checks
echo "Database: $DB";
#
if [[ "$SPLIT" == 'Sp' || "$SPLIT" == 'Sd' || "$SPLIT" == 'St' ]]; then
    echo "Split type:" $SPLIT
else
    echo "Split not specified"
    exit
fi
#
if [[ "$SUBSAMPL" == true ]]; then
    subsampling='subsampl'
    echo "Running script with subsampling"
elif [[ "$SUBSAMPL" == false ]]; then
    subsampling='nosubsampl'
    echo "Running script without subsampling"
else
    echo "No subsampling method specified"
    exit
fi


### start
echo "test: $subsampling"

DIR_DB=/DTINet/$DB

# check data folder
DIR_DATA=/DTINet/data
if [ -d "$DIR_DATA" ]; then
    echo "data folder exist, removing previous files in data"
    rm $DIR_DATA/*.txt
else
    echo "creating data folder"
    mkdir $DIR_DATA
fi

if [ -d splits ]; then
    echo "splits folder exist, removing previous files in data"
    rm splits/*.txt
else
    echo "creating splits folder"
    mkdir splits
fi


# Initialize network
echo "Removing previous files in /network and /feature"
DIR_NETWORK=/DTINet/network
[ -d "$DIR_NETWORK" ] && rm $DIR_NETWORK/*.txt
DIR_FEAT=/DTINet/feature
[ -d "$DIR_FEAT" ] && rm $DIR_FEAT/*.txt

# Moving Data
echo "Copy data from $DIR_DB"
# as normal method
cp $DIR_DB/mat*.txt $DIR_DATA/
cp $DIR_DB/Sim*.txt $DIR_DATA/
# adding splits to split folder
folder=Index_"$SPLIT"_"$subsampling"
echo "path folder:" $DIR_DB/$folder
echo  "Coping splits from"  "$DIR_DB/$folder/*.txt"
cp $DIR_DB/$folder/*.txt splits/



echo "Launching DTINet..."
echo "Compute similarity networks"
matlab -nodisplay -nosplash -nodesktop -r "run('src/compute_similarity.m');exit;"

echo "Run DCA"
matlab -nodisplay -nosplash -nodesktop -r "run('src/run_DCA.m');exit;"

# modify the log out
out_file=log_DTINet_"$DB"_"$SPLIT"_"$subsampling.out"

echo "Launching model with nohup, output in $out_file"
# test
nohup matlab -nodisplay -nosplash -nodesktop -r "tic;run('src/run_DTINet.m');toc;exit;" > $out_file 2>/dev/null & 
