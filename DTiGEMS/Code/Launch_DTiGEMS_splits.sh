
SPLIT=Sp
DATA_NAME=BIOSNAP
# DATA_PATH=/home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/$DATA_NAME
DATA_PATH=/home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/$DATA_NAME
eval "DOCKER_ID=$( docker run -d -t dtgems:1.0 bash)";

docker exec $DOCKER_ID pip install tqdm
docker exec $DOCKER_ID mkdir -p /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/
docker exec $DOCKER_ID mkdir -p /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker exec $DOCKER_ID mkdir -p /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/

docker cp $DATA_PATH/${DATA_NAME}_drug_adjmat_FDA_aers_bit.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_adjmat_FDA_aers_freq.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_Rchemcpp_lambda.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_Rchemcpp_marginalized.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_Rchemcpp_minmaxTanimoto.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_Rchemcpp_spectrum.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_Rchemcpp_tanimoto.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/${DATA_NAME}_drug_SIDER_SideEffect.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp $DATA_PATH/Drugs_SIMCOMP_scores.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/

docker cp $DATA_PATH/${DATA_NAME}_prot_GO_PPI.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_mismatch_kernel_k3m1.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_mismatch_kernel_k3m2.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_mismatch_kernel_k4m1.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_mismatch_kernel_k4m2.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_SmithWaterman_scores_MinMax.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_spectrum_kernel_k3.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp $DATA_PATH/${DATA_NAME}_prot_spectrum_kernel_k4.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/

docker exec $DOCKER_ID find /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/  -name "${DATA_NAME}_*" -exec  basename {} \; > $DATA_PATH/allTsim_files.txt
docker cp $DATA_PATH/allTsim_files.txt $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker exec $DOCKER_ID find /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/  -name "${DATA_NAME}_*" -exec  basename {} \; > $DATA_PATH/allDsim_files.txt
docker cp $DATA_PATH/allDsim_files.txt $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/

# python3 get_admat_from_dti_Yamanashi.py ./../../DB/Data/Yamanashi_et_al_GoldStandard/$DATA_NAME/interactions/${DATA_NAME,,}_admat_dgc_mat_2_line.txt
# python3 get_admat_from_dti_Davis.py ./../../DB/Data/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al_w_labels.tsv
docker cp $DATA_PATH/${DATA_NAME}_admat.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/admat_dgc.txt
docker cp $DATA_PATH/${DATA_NAME}_dti.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/DTI.txt

docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/Model_Sp_Sd_St_split_Improved.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Model_Sp_Sd_St_split_Improved.py

docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/load_datasets.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/load_datasets.py
docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/Clean_data_${DATA_NAME}.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Clean_data.py
# docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/Clean_data_Davis.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Clean_data.py
# docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/DTIs_Main.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/DTIs_Main.py
docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/DTIs_Main_Splits.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/DTIs_Main.py
docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/n2v_mainFunctions.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/n2v_mainFunctions.py

docker exec -w /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec $DOCKER_ID python3 -u /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Clean_data.py
nohup docker exec -w /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec $DOCKER_ID python3 -u /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/DTIs_Main.py --mode $SPLIT> /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Results/${DATA_NAME}_results.txt &

docker stop $DOCKER_ID && docker rm $DOCKER_ID
