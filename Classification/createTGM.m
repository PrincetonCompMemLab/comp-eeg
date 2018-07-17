function [krTraj, krLabels, skippedItems, ...
    responseTraj, winTime] = createTGM(sub, preproc, compWinToUse, ...
    dataRoot, behaveDataRoot, resultRoot, answerFile, isStudy)
% createTGM: creates a temporal generalization metrix (TGM) training on
% session 1 data and testing on session 2 data.
%
% Inputs:
%   sub: the subject to process, e.g. 'AA'
%   compWinToUse: the windows from session 1 to use for training (9:36 in
%   paper)
%   preproc: the preprocessing applied to the competition data (should
%   match that used to load the krData)
%   dataRoot: directory containing the data to be loaded
%   behaveDataRoot: directory containing behavioral data
%   resultRoot: directory where results should be stored
%   answerFile: file containing correct answers to KR questions
%   isStudy: If 1, use study trials instead of recall trials
%
% Outputs:
%   krTraj: items x rounds x sess2time x sess1time matrix of class
%   probabilities from training on session 1 and testing on session 2
%   krLabels: labels for which items were correct/incorrect in session 3
%   skippedItems: items that did not have 4 rounds of EEG data (due to
%   artifacts)
%   responseTraj: subject's responses for each item/round
%   winTime: timepoints (relative to stimulus onset) corresponding to
%   windows used for classification

if strcmp(preproc, 'theta')
    proc_str = '_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time';
elseif strcmp(preproc, 'voltage')
    proc_str = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time';
else
    proc_str = preproc;
end

if isStudy
    study_str = '_study';
else
    study_str = '';
end

KRanswer = importdata(answerFile);

numCompWin = length(compWinToUse);

loadFname = [dataRoot sub '/CompEEG__KR_' sub proc_str '.mat'];
fname = [resultRoot sub '/KRanalysis_TGM_' preproc study_str '.mat'];

krTraj = [];
krLabels = [];
load(loadFname);
featData = double(featData);
%% Predict KR

itemTraj = nan(length(KRanswer), 4, size(featData,2)/64, numCompWin);
for iWin = 1:numCompWin
    itemTraj(:,:,:,iWin) = runSess2Prediction(sub, featData, labels, compWinToUse(iWin), ...
        preproc, dataRoot, resultRoot, isStudy);
end

%% Combine with final quiz answers
load([resultRoot '/answers/' sub '_answers.mat']);
numQ = length(answerList);
corrAnswers = false(numQ, 1);
for a = 1:numQ
    corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
end

itemTrajCorr = itemTraj(corrAnswers, :, :, :);
itemTrajInc = itemTraj(~corrAnswers, :, :, :);

skippedItems = [];
indCorr = 1;
indInc = 1;
for i = 1:length(corrAnswers)
    if corrAnswers(i)
        if ~any(isnan(itemTrajCorr(indCorr,:,:,:)))
            krTraj = cat(1, krTraj, itemTrajCorr(indCorr,:,:,:));
            krLabels = cat(1, krLabels, 1);
        else
            skippedItems = cat(1, skippedItems, i);
        end
        indCorr = indCorr + 1;
    else
        if ~any(isnan(itemTrajInc(indInc,:,:,:)))
            krTraj = cat(1, krTraj, itemTrajInc(indInc,:,:,:));
            krLabels = cat(1, krLabels, 0);
        else
            skippedItems = cat(1, skippedItems, i);
        end
        indInc = indInc + 1;
    end
    
end
load([behaveDataRoot '/' sub '/' sub '_answerTraj.mat']);

numPotTraj = size(responseTraj, 1);
newResponseTraj = [];
responseTraj(isnan(responseTraj)) = 0; %#ok<*SAGROW>
for iPTraj = 1:numPotTraj
    if ~any(iPTraj == skippedItems);
        newResponseTraj = cat(1, newResponseTraj, responseTraj(iPTraj,:));
    end
end
responseTraj = newResponseTraj;
save(fname, 'krTraj', 'krLabels', 'skippedItems', 'responseTraj', 'winTime');

end