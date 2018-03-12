function aggPermResultsComp(subjects, ...
    windows, preproc, numPerm, numFolds, doZ, resultRoot)
% aggPermResultsComp: aggregates competition results over subjects
%
% Inputs:
%   subjects: the subjects to be processed, e.g. {'AA', 'BB', ...}
%   windows: the time windows to process. Note that after feature
%   extraction the data is in the form of averages over windows set during
%   the preprocessing stage. This points to which window to use as features
%   preproc: string indicating what preprocessing has been applied to the
%   data (theta or voltage)
%   numPerm: the number of features to run (100 in paper)
%   doZ: whether to zscore or not (1 in paper)
%   resultRoot: where you want the result saved. Expecting subdir for each 
%   subject

loadFileString = [resultRoot '%s/CompEEG_%s_' num2str(numFolds) ...
                    'FCV_win%d_permAccs' num2str(numPerm) ...
                    '_Z' num2str(doZ) '_' preproc '.mat'];
saveFileString = [resultRoot 'CompEEG_' num2str(numFolds) ...
                    'FCV_win%d_permAccs' num2str(numPerm) ...
                    '_Z' num2str(doZ) '_' preproc '.mat'];
%%
trueSubAccs = nan(length(subjects), numFolds);
permSubAccs = nan(length(subjects), numPerm, numFolds);
for i_w = 1:length(windows)
    w = windows(i_w);
    for s = 1:length(subjects)
        sub = subjects{s};
        loadFile = sprintf(loadFileString, sub, sub, w);
        load(loadFile);
        trueSubAccs(s, :) = trueAcc;
        permSubAccs(s, :, :) = permAccs;
    end
    saveFile = sprintf(saveFileString, w);
    save(saveFile, 'trueSubAccs', 'permSubAccs', 'winTime');
end

end