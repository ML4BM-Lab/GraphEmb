# !/bin/bash
# Launch_DTI2Vec.sh -path ./../Data/Yamanashi_et_al_GoldStandard/E/ -docker_name dti2vec:1.0

while getopts path:docker_name:f: flag
do
    case "${flag}" in
        path) path_2_data=${OPTARG};;
        docker_name) docker_name=${OPTARG};;
        f) fullname=${OPTARG};;
    esac
done
echo "path: $path_2_data";
echo "docker_name: $docker_name";
echo "Full Name: $fullname";


# This works for UXIA
# while getopts "p:d:" opt; do
#   case $opt in
#     d) dockerName=$OPTARG      ;;
#     p) path=$OPTARG   ;;
#     *) echo 'error' >&2
#        exit 1
#   esac
# done
# echo "dockerName: $dockerName";
# echo "path: $path";



# 1.Create the container from the DTI2Vec image
eval "DOCKER_ID=$( docker run -d -t $docker_name )";
echo $DOCKER_ID

# 2.Create the folders for our DB
docker exec $DOCKER_ID mkdir -p /DTi2Vec/Input/Custom/
docker exec $DOCKER_ID mkdir -p /DTi2Vec/EMBED/Custom/
docker exec $DOCKER_ID mkdir -p /DTi2Vec/test_predicted_DT/Custom/
# 3.Copy the data to the docker container
# In order to run this method, we must modify the naming of the files 
# to modify the less the original script of DTI2Vec
# The dataset will be called Custom and the needed files are:
# 	Custom_admat_dgc.txt <- adjacency matrix Drug-Prot
#	R_Custom.txt <- dti list
#  EmbeddingFold_.txt <-  10 fold embeddings on EMBED folder
docker exec $DOCKER_ID mkdir -p /DTi2Vec/Input/Custom/
docker cp ./../Data/Yamanashi_et_al_GoldStandard/E/E_admat.tsv $DOCKER_ID:/DTi2Vec/Input/Custom/Custom_admat_dgc.txt
docker cp ./../Data/Yamanashi_et_al_GoldStandard/E/E_dti.tsv $DOCKER_ID:/DTi2Vec/Input/Custom/R_Custom.txt

for ((ind = 1; ind<=10; ind++))
do 
	eval "docker cp ./../Data/Yamanashi_et_al_GoldStandard/E/E_Nodes_Embedding_64_$ind.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
done



#  docker exec -it $DOCKER_ID bash
time docker exec  $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier xgbc --func Hadmard > results_custom.txt
# docker cp  $DOCKER_ID:/DTi2Vec/results_custom.txt .


https://raw.githubusercontent.com/MahaThafar/Affinity2Vec/a242fc8e09da7f14d1d8b694a36cef7cae13a439/training_functions.py