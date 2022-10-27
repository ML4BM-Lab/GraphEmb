eval "DOCKER_ID=$( docker run -d -t dtgems:1.0 bash)";

docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/



docker exec $DOCKER_ID mkdir -p /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/
docker exec $DOCKER_ID mkdir -p /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker exec $DOCKER_ID mkdir -p /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
#
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_adjmat_FDA_aers_bit.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_adjmat_FDA_aers_freq.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_Rchemcpp_lambda.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_Rchemcpp_marginalized.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_Rchemcpp_minmaxTanimoto.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_Rchemcpp_spectrum.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_Rchemcpp_tanimoto.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_drug_SIDER_SideEffect.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/Drugs_SIMCOMP_scores.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
#
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_GO_PPI.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_mismatch_kernel_k3m1.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_mismatch_kernel_k3m2.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_mismatch_kernel_k4m1.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_mismatch_kernel_k4m2.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_SmithWaterman_scores_MinMax.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_spectrum_kernel_k3.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_prot_spectrum_kernel_k4.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
#
docker exec $DOCKER_ID ls -A1 /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/ > /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/allTsim_files.txt
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/allTsim_files.txt $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Tsim/
docker exec $DOCKER_ID ls -A1 /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/ > /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/allDsim_files.txt
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/allDsim_files.txt $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/Dsim/
#
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_admat.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/admat_dgc.txt
docker cp /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/Yamanashi_et_al_GoldStandard/E/E_dti.tsv $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/DTI.txt
#
docker cp  /home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Code/load_datasets.py $DOCKER_ID:/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/load_datasets.py

