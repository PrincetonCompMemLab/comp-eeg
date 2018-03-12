function [trueAcc, permAccs, winTime] = runCompClass_PermTest(sub, ...
    winToUse, preproc, numPerm, numFolds, doZ, dataRoot, ...
    resultRoot, randSeedFile)
% runCompClass_PermTest: a function for performing the high/low competition
% classification task for each time window, and the accompanying
% permutation test. Saves results
%
% Inputs:
%   sub: the subject to be processed, e.g. 'AA'
%   winToUse: the time window to process. Note that after feature
%   extraction the data is in the form of averages over windows set during
%   the preprocessing stage. This points to which window to use as features
%   preproc: string indicating what preprocessing has been applied to the
%   data (theta or voltage)
%   numPerm: the number of features to run (100 in paper)
%   doZ: whether to zscore or not (1 in paper)
%   dataRoot: the location of the data. Expecting subdir for each subject
%   resultRoot: where you want the result saved. Expecting subdir for each 
%   subject
%   randSeedFile: to keep the cross-validation folds consistent over
%   windows and permutations, save a random seed in this file for use
%
% Outputs:
%   trueAcc: the accuracy at this window with true label assignments
%   permAccs: the numPerm x 1 histogram of permutation accuracies
%   winTime: the time corresponding to each timepoint (relative
%       to stimulus onset).

if strcmp(preproc, 'theta')
    proc_str = '_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time';
elseif strcmp(preproc, 'voltage')
    proc_str = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time';
else
    proc_str = preproc;
end

loadFname = [dataRoot sub '/CompEEG_' sub proc_str '.mat'];
if ~exist([resultRoot sub '/'], 'dir')
    mkdir([resultRoot sub '/'])
end
resultFname = [resultRoot sub '/CompEEG_' sub '_' num2str(numFolds) ...
    'FCV_win' num2str(winToUse) '_permAccs' num2str(numPerm) ...
    '_Z' num2str(doZ) '_' preproc '.mat'];

if exist(randSeedFile, 'file')
    load(randSeedFile);
else
    randSeed = 'shuffle';
end

if ~exist(resultFname, 'file');
    if ~exist(loadFname, 'file')
        fprintf('Subject %s failed.\n', sub);
    else
        load(loadFname);
        Y = labels(:,1);
        [N, W] = size(featData);
        numFeatWins = W/64; % There are 64 channels
        
        startInd = (winToUse-1)*64 + 1;
        endInd = winToUse*64;
        X = featData(:,startInd:endInd);
        disp(size(X, 1))
        rng(randSeed);
        % Data Splitting held constant
        folds = crossvalind('Kfold', N, numFolds);
%% True Accuracy
        trueAcc = nan(numFolds, 1);
        for f = 1:numFolds
            testInds = (folds == f);
            % Zscore train data and use mean/var to adjust test data
            if doZ
                [trainData, mu, sigma] = zscore(X(~testInds, :));
                testData = (X(testInds,:) - repmat(mu, sum(testInds), 1))./repmat(sigma,sum(testInds), 1);
            else
                trainData = X(~testInds, :); %#ok<*UNRCH>
                testData = X(testInds,:);
            end
            % Use warm starts if possible
            if f > 1
                [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false, B);
            else
                [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false);
            end
            P = [testData ones(sum(testInds), 1)]*B;
            P = exp(P)./(1 + exp(P));
            Yhat = double(P > 0.5);
            
            trueAcc(f) = sum(Yhat == Y(testInds))/length(Yhat);
        end
    end
%% Permutation Accuracy
    permAccs = nan(numPerm, numFolds);
    tic
    for p = 1:numPerm
        % Permute Labels
        rng('shuffle');
        Y = Y(randperm(N));
        
        for f = 1:numFolds
            testInds = (folds == f);
            % Zscore train data and use mean/var to adjust test data
            if doZ
                [trainData, mu, sigma] = zscore(X(~testInds, :));
                testData = (X(testInds,:) - repmat(mu, sum(testInds), 1))./repmat(sigma,sum(testInds), 1);
            else
                trainData = X(~testInds, :); %#ok<*UNRCH>
                testData = X(testInds,:);
            end
            % Use warm starts if possible
            if f > 1
                [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false, B);
            else
                [B, lambda] = logReg(trainData, Y(~testInds), [1 2], false);
            end
            P = [testData ones(sum(testInds), 1)]*B;
            P = exp(P)./(1 + exp(P));
            Yhat = double(P > 0.5);
            
            permAccs(p, f) = sum(Yhat == Y(testInds))/length(Yhat);
        end
    end
    save(resultFname,...
        'trueAcc', 'permAccs', 'winTime');
    fprintf('Subject %s succeeded.\n', sub);
    toc
end

end