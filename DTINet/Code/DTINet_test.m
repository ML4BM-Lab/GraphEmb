function [roc_avg, pr_avg] = DTINet(seed, nFold, interaction, drug_feat, prot_feat, dim_imc)
	% later change seed here
	% auroc lits
	AUROC = zeros(nFold, 1);
	AUPRC = zeros(nFold, 1);

	for foldID = 1 : nFold
		train_posIdx_ = load(append('../splits/train_pos_', int2str(seed),'_', int2str(foldID), '.txt'));
		train_posIdx = sub2ind(size(interaction), train_posIdx_(:,1), train_posIdx_(:,2));
		fprintf('train pos');
		train_negIdx_ = load(append('../splits/train_neg_', int2str(seed),'_', int2str(foldID), '.txt'));
		train_negIdx = sub2ind(size(interaction), train_negIdx_(:,1), train_negIdx_(:,2))
		fprintf('train neg');
		train_idx = [train_posIdx; train_negIdx];
		Ytrain = [ones(length(train_posIdx), 1); zeros(length(train_negIdx), 1)];
		fprintf('Train data: %d positives, %d negatives\n', sum(Ytrain == 1), sum(Ytrain == 0));
		fprintf('test pos');
		test_posIdx_ = load(append('../splits/test_pos_', int2str(seed),'_', int2str(foldID), '.txt'));
		test_posIdx = sub2ind(size(interaction), test_posIdx_(:,1), test_posIdx_(:,2));
		fprintf('test negs');
		test_negIdx_ = load(append('../splits/test_neg_', int2str(seed),'_', int2str(foldID), '.txt'));
		test_negIdx = sub2ind(size(interaction), test_negIdx_(:,1), test_posIdx_(:,2));

		test_idx = [test_posIdx; test_negIdx];
		Ytest = [ones(length(test_posIdx), 1); zeros(length(test_negIdx), 1)];		
		fprintf('Test data: %d positives, %d negatives\n', sum(Ytest == 1), sum(Ytest == 0));

		[I, J] = ind2sub(size(interaction), train_idx);
		Xtrain = sparse(I, J, Ytrain, size(interaction, 1), size(interaction, 2));

		[W, H, ~] = train_mf(Xtrain, sparse(drug_feat), sparse(prot_feat), ...
						[' -l ' num2str(1) ' -k ' num2str(dim_imc) ' -t 10' ' -s ' num2str(10)]); 
		Zscore = drug_feat * W' * H * prot_feat';
		Ypred = Zscore(test_idx);

		[trainroc, trainpr] = auc(Ytrain, Zscore(train_idx), 1e-6);
		[testroc, testpr] = auc(Ytest, Ypred, 1e-6);
		AUROC(foldID) = testroc;
		AUPRC(foldID) = testpr;
		fprintf('Fold %d, Train: AUROC=%f AUPR=%f; Test: AUROC=%f, AUPR=%f\n', foldID, trainroc, trainpr, testroc, testpr);
	end
	roc_avg = mean(AUROC);
    pr_avg = mean(AUPRC);
end
