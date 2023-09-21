# 6 DRUG/PROTEIN RELATED MATRICES

# DRUG-PROTEIN (1) -> DRUGBANK
# DRUG-DRUG (1) -> DRUGBANK
# PROTEIN-PROTEIN (1) -> HPRD DATABASE (UXIA)
# DRUG-DISEASE / PROTEIN-DISEASE (2) -> COMPARATIVE TOXICOGENOMICS DATABASE (UXIA)
# DRUG-SIDE-EFFECT (1) -> SIDER ✓

## --- ##
# Drug Chemical Structure Information -> Pairwise matrix built with the Morgan Fingerprints with radius 2 (RDKIT) ✓✓
# Protein Sequence Information -> Pairwise Smith-Waterman ✓

## --- HOW TO RUN ---- ##

# -d: The embedding dimension d, default: 1024.
# -n: Global norm to be clipped, default: 1.
# -k: The dimension of project matrices, default: 512.
# -r: Positive and negative. Two choices: ten and all, the former one sets the positive:negative = 1:10, the latter one considers all unknown DTIs as negative examples. Default: ten.
# -t: Test scenario. The DTI matrix to be tested. Choices are: o, mat_drug_protein.txt will be tested; homo, mat_drug_protein_homo_protein_drug.txt will be tested; drug, mat_drug_protein_drug.txt will be tested; disease, mat_drug_protein_disease.txt will be tested; sideeffect, mat_drug_protein_sideeffect.txt will be tested; unique, mat_drug_protein_drug_unique.txt will be tested. Default: o.

# This is the modified NeoDTI
nohup python -u src/NeoDTI_cv_mod.py -m Yamanishi_et_al_GoldStandard/NR > logs/NR_default.out &
nohup python -u src/NeoDTI_cv_mod.py -m Yamanishi_et_al_GoldStandard/IC > logs/IC_splits_Sp.out &

nohup python -u src/NeoDTI_cv_mod.py -m Yamanishi_et_al_GoldStandard/GPCR > logs/GPCR_splits_St.out &

docker cp nr_splits.pickle sad_hodgkin:/NeoDTI/data/Yamanishi_et_al_GoldStandard/NR/.


# -- EXTRA
dic_smiles.json
disease.txt
drug.txt
DTI_*.txt
protein.txt

# -- DRUG MATRICES
## mat_drug_drug (DDI)
## Similarity_Matrix_Drugs (Drug Properties)
## mat_drug_disease (diseases)
## mat_drug_se (side effects)

# -- DTI matrices
## mat_drug_protein.txt (dti)

# -- PROTEIN MATRICES
## mat_protein_protein (PPI)
## Similarity_Matrix_Proteins (Protein Properties)
## mat_protein_disease # (diseases)






