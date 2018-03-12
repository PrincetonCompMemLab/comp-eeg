function windowMaps = weight2map(subjects, window, doZ, preproc, dataRoot, resultRoot)
% weight2map: converts logistic regression weights to interpretable map
% based on the technique from Haufe et al. 2014. Plots the given window as
% a topoplot
%
% Inputs:
%   subjects: subjects to process (e.g. {'AA', 'BB', ...}
%   window: window of timepoints to average maps over for topoplot
%   preproc: string indicating what preprocessing has been applied to the
%   data (theta or voltage)
%   doZ: whether to zscore or not (1 in paper)
%   dataRoot: the location of the data. Expecting subdir for each subject
%   resultRoot: where you want the result saved. Expecting subdir for each 
%   subject
% Outputs:
%   windowMaps: map for the requested window

loc_file = '/Users/nrafidi/Documents/MATLAB/compEEG-data/electrode_locs.xyz';

if strcmp(preproc, 'theta')
    proc_str = '_Vis_BP2-200_N60_Ref_Hilbert-theta_Epochs_Features_Overlap_Time';
elseif strcmp(preproc, 'voltage')
    proc_str = '_Vis_BP2-200_N60_Ref_Epochs_Base_ICA1-2_Features_Overlap_Time';
else
    proc_str = preproc;
end
%%
load([resultRoot 'CompEEG_Models_Z' num2str(doZ) ...
    '_' preproc '.mat']);

[num_sub, ~] = size(subModels);
num_win = size(subModels{1}, 1);
%%
mapFile = [resultRoot 'haufeMaps_Z' num2str(doZ) '_' preproc '.mat'];

if ~exist(mapFile, 'file')
    subMaps = nan(num_sub, num_win, 65);
    for i_sub = 1:num_sub
        sub = subjects{i_sub};
        subW = subModels{i_sub};
        loadFname = [dataRoot sub '/CompEEG_' sub proc_str '.mat'];
        load(loadFname);
        for p = 1:num_win

            W = subW(p, :)';

            startInd = (p-1)*64 + 1;
            endInd = p*64;
            X = [zscore(featData(:,startInd:endInd)), ones(size(featData, 1), 1)];

            S = X*W;

            cov_X = cov(X);

            cov_S = cov(S);

            subMaps(i_sub, p, :) = cov_X*W/cov_S;

        end
    end

    save(mapFile, 'subMaps');
else
    load(mapFile);
end

windowMaps = squeeze(mean(mean(subMaps(:, window, :), 1), 2));

absMax = max(abs(windowMaps));


eeglab;

map = [windowMaps(1:64); zeros(8, 1)];
f = figure;
topoplot(map, loc_file, 'maplimits', [-absMax, absMax]);
colorbar
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f, [resultRoot 'figures/CompEEG_haufeTopo_' num2str(min(window)) ...
    '-' num2str(max(window)) '_Z' num2str(doZ) '_' preproc '.pdf']);

end