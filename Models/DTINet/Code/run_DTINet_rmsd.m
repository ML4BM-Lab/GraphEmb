dim_drug = 100;
dim_prot = 400;
dim_imc = 50;
rm_homo = false;
if rm_homo == true;
	interaction = load('../data/mat_drug_protein_remove_homo.txt');
else
	interaction = load('../data/mat_drug_protein.txt');
end
drug_feat = load(['../feature/drug_vector_d',num2str(dim_drug),'.txt']);
prot_feat = load(['../feature/protein_vector_d',num2str(dim_prot),'.txt']);
nFold = 10;
Nrepeat = 1; % only 1 for rmsd option now

AUROC_test = zeros(Nrepeat, 1);
AUPRC_test = zeros(Nrepeat, 1);

AUROC_final = zeros(Nrepeat, 1);
AUPR_final = zeros(Nrepeat, 1);

for p = 1 : Nrepeat
    fprintf('Repetition #%d\n', p);
    [AUROC_test(p), AUPRC_test(p), AUROC_final(p), AUPR_final(p)] = DTINet(p, nFold, interaction, drug_feat, prot_feat, dim_imc);
end

for i = 1 : Nrepeat
	fprintf('Repetition #%d: AUROC=%.6f, AUPR=%.6f\n', i, AUROC_test(i), AUPRC_test(i));
end

fprintf('TEST. Mean: AUROC=%.6f, AUPR=%.6f\n', mean(AUROC_test), mean(AUPRC_test));
fprintf('Final Fold. Mean: AUROC=%.6f, AUPR=%.6f\n', mean(AUROC_final), mean(AUPR_final));