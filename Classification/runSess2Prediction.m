function [ itemTraj ] = runSess2Prediction(sub, ...
    krData, krLabels, compWin, preproc, dataRoot, resultRoot, isStudy)
% runSess2Prediction: trains a classifier on session 1 data from the
% given window and applies it to all windows in the given session 2 data
%
% Inputs:
%   sub: the subject to process, e.g. 'AA'
%   krData: the subject's session 2 data
%   krLabels: the labels for the subject's session 2 data
%   compWin: the window in the session 1 experiment to train on
%   preproc: the preprocessing applied to the competition data (should
%   match that used to load the krData)
%   dataRoot: directory containing the data to be loaded
%   resultRoot: directory where results should be stored
%   isStudy: If 1, use study trials instead of recall trials
%
% Outputs:
%   itemTraj: items x rounds x time matrix containing the class probability
%   on the session 2 data for each item, round, and timepoint
numChan = 64;

if isStudy
    label_rows = 0;
else
    label_rows = 1;
end

if strcmp(preproc, 'theta')
    proc_str = '_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time';
elseif strcmp(preproc, 'voltage')
    proc_str = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time';
else
    proc_str = preproc;
end

fnameWeights = [resultRoot sub '/CompEEG_Weights_cWin' num2str(compWin) '_' preproc '.mat'];

rng('shuffle');

compFname = [dataRoot sub '/CompEEG_' sub proc_str '.mat'];

if exist(fnameWeights, 'file')
    load(fnameWeights);
else
    load(compFname);
    
    featData = double(featData);
    [compData, mu, sigma] = zscore(featData);
    compLabels = labels(:,1); %#ok<*NODEF>
    
    compWinToUse = (compWin*numChan + 1):((compWin+1)*numChan);
    
    B = logReg(compData(:,compWinToUse), compLabels, [1 10], false);
    
    save(fnameWeights, 'B', 'mu', 'sigma');
end

krData = krData(krLabels(:,1) == label_rows, :);
krLabels = krLabels(krLabels(:,1) == label_rows, 2);


numWinToTry = size(krData,2)/numChan;
uniqueItems = unique(krLabels);
if any(uniqueItems == 255)
    uniqueItems = uniqueItems(uniqueItems ~= 255);
end
itemTraj = nan(length(uniqueItems), 4, numWinToTry);
for i = 1:length(uniqueItems)
    
    items = krLabels == (i+1);
    for w = 1:numWinToTry
        timeWindow = ((w-1)*numChan + 1):(w*numChan);
        krData(items, timeWindow) = (krData(items,timeWindow) - repmat(mu(timeWindow), sum(items), 1))./repmat(sigma(timeWindow), sum(items), 1);
        P = [krData(items,timeWindow) ones(sum(items), 1)]*B;
        probs = exp(P)./(1+exp(P));
        if length(probs) == 4
            itemTraj(i, :, w) = probs';
        end
    end
    
end

end

