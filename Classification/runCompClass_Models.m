function [models, winTime] = runCompClass_Models(sub, preproc, doZ, ...
    dataRoot, resultRoot)
% runCompClass_Models: a function that trains a logistic regression model
% to distinguish high/low competition states on all the data and saves the
% learned weights, also returning them.
%
% Inputs:
%   sub: the subject to be processed, e.g. 'AA'
%   preproc: string indicating what preprocessing has been applied to the
%   data (theta or voltage)
%   doZ: whether to zscore or not (1 in paper)
%   dataRoot: the location of the data. Expecting subdir for each subject
%   resultRoot: where you want the result saved. Expecting subdir for each 
%   subject
%
% Outputs:
%   models: one model for each timepoint
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

if ~exist(loadFname, 'file')
    fprintf('Subject %s failed.\n', sub);
else
    load(loadFname);
    
    Y = labels(:,1);
    [~, W] = size(featData);
    numFeatWins = W/64;
    models = cell(numFeatWins, 1);
    
    for p = 1:numFeatWins
        startInd = (p-1)*64 + 1;
        endInd = p*64;
        X = featData(:,startInd:endInd);
        
        if doZ
            [trainData, ~, ~] = zscore(X);
        else
            trainData = X;
        end
        
        [B, ~] = logReg(trainData, Y, [1 2], false);
        
        
        models{p} = B;
        
    end
    save([resultRoot sub '/CompEEG_' sub '_Models_Z' num2str(doZ) ...
        '_' preproc '.mat'], 'models', 'winTime');
end

end