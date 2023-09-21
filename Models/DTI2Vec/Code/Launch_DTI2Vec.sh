# !/bin/bash
# Launch_DTI2Vec.sh -p ./../Data/Yamanashi_et_al_GoldStandard/NR -d dti2vec:1.0 
# Launch_DTI2Vec.sh -p ./../Data/BIOSNAP -d dti2vec:1.0 
# usage() { echo "Usage: $0 [-s <45|90>] [-p <string>]" 1>&2; exit 1; }

# nohup bash Launch_DTI2Vec.sh -p ./../Data/BIOSNAP -d dti2vec:1.0 &
# nohup bash Launch_DTI2Vec.sh -p ./../Data/DrugBank -d dti2vec:1.0 &
# nohup bash Launch_DTI2Vec.sh -p ./../Data/BindingDB -d dti2vec:1.0 &
# nohup bash Launch_DTI2Vec.sh -p ./../Data/Davis_et_al -d dti2vec:1.0 -m Sp &
# nohup bash Launch_DTI2Vec.sh -p ./../Data/Yamanashi_et_al_GoldStandard/NR -d dti2vec:1.0 -m Sp &
# nohup bash Launch_DTI2Vec.sh -p ./../Data/BindingDB -d dti2vec:1.0 -m Sp &

mode='DEFAULT'
while getopts "p:d:m:" opt; do
    echo $opt
    case $opt in
        p) path_2_data=$OPTARG;;
        d) docker_name=$OPTARG;;
        m) mode=$OPTARG;;
        *) echo 'Error in command line parsing' >&3
           exit 1
    esac
done
echo "mode: $mode"
echo "path: $path_2_data";
echo "docker_name: $docker_name";


#====================== DOCKER PREPARATION======================
# 1.Create the container from the DTI2Vec image
eval "DOCKER_ID=$( docker run -d -t $docker_name bash)";
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
    GS=false
    mkdir -p ./../Results/Yamanashi_et_al_GoldStandard/$sub_db
    echo "Results would be saved in: ./../Results/Yamanashi_et_al_GoldStandard/$sub_db"
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



if [ $GS = true ]; then
    echo 'Non Yamanishi-like'
    echo "GS: $GS"
    if [[ $path_2_data =~ (Data\/)([A-Za-z_]*) ]]; then
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
    # k is the number of 'augmented' neighbors added
    for k in 2 5 
    do echo "k: $k";
        # embeddings is the number of dimensions created for each node in N2Vec
        for embeddings in 128_ 256_ 512_ 
        do echo ""embeddings: $embeddings"";
            for ((ind = 1; ind<=10; ind++))
            do 
                echo "docker cp $path_2_data/$DB$emb_file$embeddings${k}_$ind.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
                eval "docker cp $path_2_data/$DB$emb_file$embeddings${k}_$ind.emb $DOCKER_ID:/DTi2Vec/EMBED/Custom/EmbeddingFold_$ind.txt"
            done
            for func in Hadmard Concat
            do
                if [ $mode == "DEFAULT" ]; then
                        mkdir -p ./../Results/$DB/
                        echo "Results would be saved in: ./../Results/$DB/results_custom_${k}_$embeddings$func.txt"
                        echo "==========================================================="
                        echo "class: xgbc";
                        echo "func: $func";
                        echo "Grid Search: $GS";
                        echo "==========================================================="
                        echo "docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier xgbc --func $func > ./../Results/$DB/results_custom_${k}_$embeddings$func.txt"
                        docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier xgbc --func $func > ./../Results/$DB/results_custom_${k}_$embeddings$func.txt
                fi
                if [ $mode != 'DEFAULT' ]; then
                        mkdir -p ./../Results/$DB/
                        docker exec  $DOCKER_ID pip3 install tqdm
                        # modified utilities script to parse the split modes
                        echo "docker cp ./load_datasets_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
                        eval "docker cp ./load_datasets_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
                        # main modified to create the splits
                        echo "docker cp ./DTi2vec_main_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
                        eval "docker cp ./DTi2vec_main_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
                        # the script where the splits are defined
                        echo "docker cp ./Model_Sp_Sd_St_split_Improved.py $DOCKER_ID:/DTi2Vec"
                        eval "docker cp ./Model_Sp_Sd_St_split_Improved.py $DOCKER_ID:/DTi2Vec"
                        echo "Results would be saved in: ./../Results/$DB/results_custom_${k}_$embeddings$func.txt"
                        echo "==========================================================="
                        echo "class: xgbc";
                        echo "func: $func";
                        echo "Grid Search: $GS";
                        echo "Split mode: $mode";
                        echo "==========================================================="
                        echo "docker exec $DOCKER_ID python3 -u ./DTi2vec_main_SD_ST_SP.py --data Custom --classifier xgbc --func $func --mode $mode > ./../Results/$DB/results_custom_${k}_$embeddings$func.txt"
                        docker exec $DOCKER_ID python3 -u ./DTi2vec_main_SD_ST_SP.py --data Custom --classifier xgbc --func $func --mode $mode > ./../Results/$DB/results_custom_SP_SD_ST_${k}_$embeddings$func$mode.txt
                fi
            done
        done
    done
else
    # Yamanashi datasets without Grid search
    echo "Yamanishi like"
    if [ $mode == "DEFAULT" ]; then
        echo "docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier $class --func $func > /DTi2Vec/results_custom.txt"
        docker exec $DOCKER_ID python3 -u ./DTi2vec_main.py --data Custom --classifier $class --func $func > ./../Results/Yamanashi_et_al_GoldStandard/$sub_db/results_custom.txt
    fi
    if [ $mode != 'DEFAULT' ]; then
        docker exec  $DOCKER_ID pip3 install tqdm
        echo "docker cp ./load_datasets_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
        eval "docker cp ./load_datasets_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
        # main modified to create the splits
        echo "docker cp ./DTi2vec_main_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
        eval "docker cp ./DTi2vec_main_SD_ST_SP.py $DOCKER_ID:/DTi2Vec"
        # the script where the splits are defined
        echo "docker cp ./Model_Sp_Sd_St_split_Improved.py $DOCKER_ID:/DTi2Vec"
        eval "docker cp ./Model_Sp_Sd_St_split_Improved.py $DOCKER_ID:/DTi2Vec"
        echo "==========================================================="
        echo "class: $class";
        echo "func: $func";
        echo "Grid Search: $GS";
        echo "Split mode: $mode";
        echo "==========================================================="
        echo "docker exec $DOCKER_ID python3 -u ./DTi2vec_main_SD_ST_SP.py --data Custom --classifier $class --func $func --mode $mode >  ./../Results/Yamanashi_et_al_GoldStandard/$sub_db/results_custom$mode.txt"
        docker exec $DOCKER_ID python3 -u ./DTi2vec_main_SD_ST_SP.py --data Custom --classifier $class --func $func --mode $mode >  ./../Results/Yamanashi_et_al_GoldStandard/$sub_db/results_custom$mode.txt
    fi
fi

echo 'Done!';
docker stop $DOCKER_ID && docker rm $DOCKER_ID

# https://raw.githubusercontent.com/MahaThafar/Affinity2Vec/a242fc8e09da7f14d1d8b694a36cef7cae13a439/training_functions.py


