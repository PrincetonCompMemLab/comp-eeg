function plotTGM_subZ(subjects, preproc, compWinToUse, ...
    resultRoot, figTag, behav_pattern)
% plotTGM_subZ: plots the temporal generalization metrix (TGM) training on
% session 1 data and testing on session 2 data. Zscores the score for each
% subject.
%
% Inputs:
%   subjects: the subjects to process, e.g. {'AA', 'BB', ...}
%   compWinToUse: the windows from session 1 to use for training (9:36 in
%   paper)
%   preproc: the preprocessing applied to the competition data (should
%   match that used to load the krData)
%   resultRoot: directory where results should be stored
%   figTag: additional naming information for plots (e.g. '_test')
%   behav_pattern: either 'all', or a matrix describing which session 2
%   behavioral patterns to plot, e.g. [0, 1, 1, 1] or [0, 0, 1, 1;0, 1, 1,
%   1]
if strcmp(behav_pattern, 'all')
    behav_str = behav_pattern;
    behav_pattern = [0, 0, 1, 1; 0, 1, 0, 1; ...
    0, 1, 1, 0; 0, 1, 1, 1; 1, 0, 0, 1; ...
    1, 0, 1, 0; 1, 0, 1, 1; 1, 1, 0, 0; ...
    1, 1, 0, 1; 1, 1, 1, 0; 1, 1, 1, 1];
elseif size(behav_pattern, 1) > 1
    behav_str = num2str(reshape(behav_pattern', 1, []));
else
    behav_str = num2str(behav_pattern);
end
behav_str(isspace(behav_str)) = [];

numSub = length(subjects);

fnameString = [resultRoot '%s/KRanalysis_TGM_' preproc '.mat'];
figureString = [resultRoot 'figures/KR_TGM_FirstCorrect-R4_%s_' preproc '_' ...
    behav_str '_subZ' figTag '.pdf'];
%%
diffMatR = [];
diffMatF = [];
firstRoundCorr_R = [];
firstRoundCorr_F = [];
IndividualSubjectFirstCorrR = cell(numSub, 1);
IndividualSubjectFirstCorrF = cell(numSub, 1);
for i = 1:numSub
    load(sprintf(fnameString, subjects{i}));
    
    goodItems = false(size(responseTraj, 1), 1);
    for i_patt = 1:size(behav_pattern, 1)
        goodItems = goodItems | ismember(responseTraj, behav_pattern(i_patt, :), 'rows');
    end
    
    for j = 1:size(responseTraj, 1)
        if goodItems(j)
            if krLabels(j) == 1
                firstRoundCorr_R = cat(1, firstRoundCorr_R, find(responseTraj(j,:), 1, 'first'));
                IndividualSubjectFirstCorrR{i} = cat(1, IndividualSubjectFirstCorrR{i}, find(responseTraj(j,:), 1, 'first'));
            else
                firstRoundCorr_F = cat(1, firstRoundCorr_F, find(responseTraj(j,:), 1, 'first'));
                IndividualSubjectFirstCorrF{i} = cat(1, IndividualSubjectFirstCorrF{i}, find(responseTraj(j,:), 1, 'first'));
            end
        end
    end
    
    krTrajToAdd_R = krTraj(krLabels == 1 & goodItems,:,:,:);
    krTrajToAdd_F = krTraj(krLabels == 0 & goodItems,:,:,:);
    [num_samp_R, ~, numK, numC] = size(krTrajToAdd_R);
    num_samp_F = size(krTrajToAdd_F, 1);
    
    krScoreToAdd_R = nan(num_samp_R, numK, numC);
    krScoreToAdd_F = nan(num_samp_F, numK, numC);
    
    for iSamp = 1:num_samp_R
        krScoreToAdd_R(iSamp,:,:) = squeeze(krTrajToAdd_R(iSamp,IndividualSubjectFirstCorrR{i}(iSamp),:,:) - ...
            krTrajToAdd_R(iSamp,4,:,:));
    end
    for iSamp = 1:num_samp_F
        krScoreToAdd_F(iSamp,:,:) = squeeze(krTrajToAdd_F(iSamp,IndividualSubjectFirstCorrF{i}(iSamp),:,:) - ...
            krTrajToAdd_F(iSamp,4,:,:));
    end
    
    krScoreTotal = zscore([krScoreToAdd_R; krScoreToAdd_F], 1);
    
    
    diffMatR = cat(1, diffMatR, krScoreTotal(1:num_samp_R, :, :));
    diffMatF = cat(1, diffMatF, krScoreTotal((num_samp_R + 1):end, :, :));
end

ctWinTime = -100:20:960;
ctWinTime = ctWinTime(compWinToUse);
krWinToUse = 6:54;
krWinTime = winTime(krWinToUse);
diffMatF = diffMatF(:,krWinToUse, :);
diffMatR = diffMatR(:,krWinToUse, :);
%% Plot desired value

numR = size(diffMatR, 1);
numF = size(diffMatF, 1);
numI = min([numR, numF]);
diffMatR = diffMatR(1:numI, :, :);
diffMatF = diffMatF(1:numI, :, :);


[~, rep_corr] = ttest(diffMatR, diffMatF, 'Tail', 'right');
rep_corr = squeeze(rep_corr);

[~, rep_opp] = ttest(diffMatR, diffMatF, 'Tail', 'left');
rep_opp = squeeze(rep_opp);

rep_diff = squeeze(mean(diffMatR - diffMatF, 1));

rep = zeros(size(rep_opp));
rep(rep_opp < 0.05) = -1;
rep(rep_opp < 0.01) = -2;
rep(rep_opp < 0.001) = -5;

rep(rep_corr < 0.05) = 1;
rep(rep_corr < 0.01) = 2;
rep(rep_corr < 0.001) = 5;
%%
f1 = figure;
imagesc(krWinTime, ctWinTime, rep_diff');
ylabel('Competition Training Time (ms)');
xlabel('Swahili Learning Test Time (ms)');
colorbar
title(sprintf('Average Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);
set(f1, 'Color', 'w');
export_fig(f1, sprintf(figureString, 'abs'));

f2 = figure;
imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
ylabel('Competition Training Time (ms)');
xlabel('Swahili Learning Test Time (ms)');
colorbar
title(sprintf('P value Difference\n%s Trials', behav_str));
set(gca, 'FontSize', 18);
set(f2, 'Color', 'w');
export_fig(f2, sprintf(figureString, 'PVal'));

end