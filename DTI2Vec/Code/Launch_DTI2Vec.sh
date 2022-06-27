# !/bin/bash
# Launch_DTI2Vec.sh -p ./../Data/Yamanashi_et_al_GoldStandard/NR -d dti2vec:1.0 
# Launch_DTI2Vec.sh -p ./../Data/BIOSNAP -d dti2vec:1.0 
# usage() { echo "Usage: $0 [-s <45|90>] [-p <string>]" 1>&2; exit 1; }

# nohup bash Launch_DTI2Vec.sh -p ./../Data/BIOSNAP -d dti2vec:1.0 &

while getopts "p:d:" opt; do
    case $opt in
        p) path_2_data=$OPTARG;;
        d) docker_name=$OPTARG;;
         *) echo 'error' >&2
       exit 1
    esac
done
echo "path: $path_2_data";
echo "docker_name: $docker_name";




#====================== DOCKER PREPARATION======================
# 1.Create the container from the DTI2Vec image
eval "DOCKER_ID=$( docker run -d -t $docker_name )";
echo "The id of the running container is: $DOCKER_ID"

# 2.Create the folders for our DB
docker exec $DOCKER_ID mkdir -p /DTi2Vec/Input/Custom/
docker exec $DOCKER_ID mkdir -p /DTi2Vec/EMBED/Custom/
docker exec $DOCKER_ID mkdir -p /DTi2Vec/test_predicted_DT/Custom/
docker exec $DOCKER_ID mkdir -p /DTi2Vec/Input/Custom/



#====================== FILES PREPARATION======================
# 3.Copy the data to the docker container
# In order to run this method, we must modify the naming of the files 
# to modify the less the original script of DTI2Vec
# The dataset will be called Custom and the needed files are:
# 	Custom_admat_dgc.txt <- adjacency matrix Drug-Prot
#	R_Custom.txt <- dti list
#  EmbeddingFold_.txt <-  10 fold embeddings on EMBED folder

GS=true
# YAMANISHI CASES
if [[ "$path_2_data" == *"Yamanashi_et_al_GoldStandard"* ]]; then
    echo "Yamanashi_et_al_GoldStandard";
    if [[ $path_2_data =~ (GoldStandard\/)([A-Z]*) ]]; then
        sub_db=${BASH_REMATCH[2]}
        echo "Yamanishi :: $sub_db";
    else
        echo 'No match';
    fi
    case $sub_db in 
        NR)
            class="ab"
            func="WL1"
            GS=false
            ;;
        GPCR)
            class="xgbc"
            func="Hadmard"
            GS=false
            ;;
        IC)
            class="xgbc"
            func="Concat"
            GS=false
            ;;
        E)
            class="xgbc"
            func="Concat"
            GS=false
            ;;
        *) echo 'No default DB used, using grid search'
            class=""
            func=""
            GS=true
            ;;
    esac
    mkdir -p ./../Results/Yamanashi_et_al_GoldStandard/$sub_db
    echo "Results would be saved in: ./../Results/Yamanashi_et_al_GoldStandard/$sub_db"
    echo "==========================================================="
    echo "class: $class";
    echo "func: $func";
    echo "Grid Search: $GS";
    echo "==========================================================="
    admat="_admat.tsv"
    dti="_dti.tsv"
    emb_file="_Nodes_Embedding_"
    embeddings="128_"
    k="2_"
    docker cp  $path_2_data/$sub_db$admat $DOCKER_ID:/DTi2Vec/Input/Custom/Custom_admat_dgc.txt
    docker cp  $path_2_data/$sub_db$dti $DOCKER_ID:/DTi2Vec/Input/Custom/R_Custom.txt
    for ((ind = 1; ind<=10; ind++))
    do 
        # echo "docker cp $path_2_data/$sub_db$emb_file$embeddings$k$ind.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
        eval "docker cp $path_2_data/$sub_db$emb_file$embeddings$k$ind.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
    done
fi


# OTHER CASES with GS



echo "here";
echo "GS: $GS"
if $GS; then
    if [[ $path_2_data =~ (Data\/)([A-Z]*) ]]; then
        DB=${BASH_REMATCH[2]}
        echo "DB :: $DB";
        else
            echo 'No match';
            exit 1;
    fi
    admat="_admat.tsv"
    dti="_dti.tsv"
    emb_file="_Nodes_Embedding_"
    docker cp  $path_2_data/$DB$admat $DOCKER_ID:/DTi2Vec/Input/Custom/Custom_admat_dgc.txt
    docker cp  $path_2_data/$DB$dti $DOCKER_ID:/DTi2Vec/Input/Custom/R_Custom.txt
    for k in 2 5 10
    do echo "k: $k";
        for embeddings in 128_ 256_ 512_ 
        do echo ""embeddings: $embeddings"";
            for ((ind = 1; ind<=10; ind++))
            do 
                echo "docker cp $path_2_data/$DB$emb_file$embeddings$k$ind_.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
                eval "docker cp $path_2_data/$DB$emb_file$embeddings$k$ind_.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
            done
            for func in Hadmard Concat
            do  mkdir -p ./../Results/$DB/
                echo "Results would be saved in: ./../Results/$DB/results_custom_${k}_$embeddings$func.txt"
                echo "==========================================================="
                echo "class: xgbc";
                echo "func: $func";
                echo "Grid Search: $GS";
                echo "==========================================================="
                echo "docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier xgbc --func $func > ./../Results/$DB/results_custom_${k}_$embeddings$func.txt"
                docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier xgbc --func $func > ./../Results/$DB/results_custom_${k}_$embeddings$func.txt
            done
        done
    done
else
   echo "docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier $class --func $func > /DTi2Vec/results_custom.txt"
   docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier $class --func $func > ./../Results/Yamanashi_et_al_GoldStandard/$sub_db/results_custom.txt
fi

echo 'Done!';

docker stop $DOCKER_ID && docker rm $DOCKER_ID

# https://raw.githubusercontent.com/MahaThafar/Affinity2Vec/a242fc8e09da7f14d1d8b694a36cef7cae13a439/training_functions.py


