% this is script has been modified for instroducing splits
function [roc_avg, pr_avg, final_roc_avg, final_pr_av] = DTINet(seed, nFold, interaction, drug_feat, prot_feat, dim_imc)
	% later change seed here
	% auroc lits
	AUROC_test = zeros(nFold, 1);
	AUPRC_test = zeros(nFold, 1);
	
	% Final Fold
	AUROC_final = zeros(nFold, 1);
	AUPRC_final = zeros(nFold, 1);

	% anhadir cargar el final fold que es el mismo para todos los folds
	% seed queda de momento porque va a ser solo 1
	finalfold_posIdx = load(append('../splits/finalfold_pos_', int2str(seed), '.txt'));
	finalfold_negIdx = load(append('../splits/finalfold_neg_', int2str(seed), '.txt'));
	final_idx = [finalfold_posIdx; finalfold_negIdx];
	Yfinal = [ones(length(finalfold_posIdx), 1); zeros(length(finalfold_negIdx), 1)];
	fprintf('Final Fold data: %d positives, %d negatives\n', sum(Yfinal == 1), sum(Yfinal == 0));

	for foldID = 1 : nFold
		train_posIdx = load(append('../splits/train_pos_', int2str(seed),'_', int2str(foldID), '.txt'));
		train_negIdx = load(append('../splits/train_neg_', int2str(seed),'_', int2str(foldID), '.txt'));
		train_idx = [train_posIdx; train_negIdx];
		Ytrain = [ones(length(train_posIdx), 1); zeros(length(train_negIdx), 1)];
		fprintf('Train data: %d positives, %d negatives\n', sum(Ytrain == 1), sum(Ytrain == 0));

		test_posIdx = load(append('../splits/test_pos_', int2str(seed),'_', int2str(foldID), '.txt'));
		test_negIdx = load(append('../splits/test_neg_', int2str(seed),'_', int2str(foldID), '.txt'));
		test_idx = [test_posIdx; test_negIdx];
		Ytest = [ones(length(test_posIdx), 1); zeros(length(test_negIdx), 1)];		
		fprintf('Test data: %d positives, %d negatives\n', sum(Ytest == 1), sum(Ytest == 0));

		[I, J] = ind2sub(size(interaction), train_idx);
		Xtrain = sparse(I, J, Ytrain, size(interaction, 1), size(interaction, 2));

		[W, H, ~] = train_mf(Xtrain, sparse(drug_feat), sparse(prot_feat), ...
						[' -l ' num2str(1) ' -k ' num2str(dim_imc) ' -t 10' ' -s ' num2str(10)]); 
		%
		Zscore = drug_feat * W' * H * prot_feat';
		% prediciones en Y para test
		Ypred = Zscore(test_idx);
		% calcula AUROC /AUPR 
		[trainroc, trainpr] = auc(Ytrain, Zscore(train_idx), 1e-6);
		[testroc, testpr] = auc(Ytest, Ypred, 1e-6);
		AUROC_test(foldID) = testroc;
		AUPRC_test(foldID) = testpr;

		% Final Fold
		Ypredfinal = Zscore(final_idx);
		[finalroc, finaltestpr] = auc(Yfinal, Ypredfinal, 1e-6);
		AUROC_final(foldID) = finalroc;
		AUPRC_final(foldID) = finaltestpr;
		
		% añadir aqui final fold
		fprintf('Fold %d, Train: AUROC=%f AUPR=%f; Test: AUROC=%f, AUPR=%f, FFold: AUROC=%f, AUPR=%f\n', foldID, trainroc, trainpr, testroc, testpr, finalroc, finaltestpr);
	end
	roc_avg = mean(AUROC_test);
    pr_avg = mean(AUPRC_test);
	final_roc_avg = mean(AUROC_final);
	final_pr_av = mean(AUPRC_final);
end
