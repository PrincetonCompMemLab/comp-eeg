function subModels = plotModelRDM(subjects, doZ, preproc, resultRoot)
% plotModelRDM: plots RDM of weight matrix averaged over subjects. Also
% savesthe aggregation of models over subjects.
%
% Inputs:
%   subjects: subjects to process (e.g. {'AA', 'BB', ...}
%   preproc: string indicating what preprocessing has been applied to the
%   data (theta or voltage)
%   doZ: whether to zscore or not (1 in paper)
%   resultRoot: where you want the result saved. Expecting subdir for each 
%   subject
% Outputs:
%   subModels: models for all input subjects
addpath /Users/nrafidi/Documents/MATLAB/Toolboxes/export_fig/

loadFileString = [resultRoot '%s/CompEEG_%s_Models_Z' num2str(doZ) ...
        '_' preproc '.mat'];
    
saveFileString = [resultRoot 'CompEEG_Models_Z' num2str(doZ) ...
    '_' preproc '.mat'];


subModels = cell(length(subjects), 1);
for s = 1:length(subjects)
    sub = subjects{s};
    loadFile = sprintf(loadFileString, sub, sub);
    load(loadFile);
    subModels{s} = cell2mat(models(:, 1)')';
end
save(saveFileString, 'subModels');

numSub = size(subModels, 1);
avgRDM = zeros(length(models));
timeVec = winTime(1:length(models));

for s = 1:numSub
    subModel = subModels{s};
    subDist = pdist(subModel, 'cosine');
    subDist = squareform(subDist);
    subDist = subDist./max(max(subDist));
    avgRDM = avgRDM + subDist;
end
avgRDM = avgRDM./numSub;
f = figure;
imagesc(timeVec, timeVec, avgRDM);
title(sprintf('Average RDM of Weight Vectors over time\n%s', preproc));
xlabel('Time Relative to Onset (ms)');
ylabel('Time Relative to Onset (ms)');
colorbar

set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f, [resultRoot 'figures/CompEEG_ModelRDM_Z' num2str(doZ) ...
    '_' preproc '.pdf']);