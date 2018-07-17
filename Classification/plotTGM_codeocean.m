function plotTGM_codeocean(subjects, preproc, compWinToUse, krWinToUse, ...
    resultRoot, figTag, clusterOptions, behav_pattern, doSubZ, isStudy)
% plotTGM: plots the temporal generalization metrix (TGM) training on
% session 1 data and testing on session 2 data. Performs cluster
% permutation test and subject-level bootstrap.
%
% Inputs:
%   subjects: the subjects to process, e.g. {'AA', 'BB', ...}
%   compWinToUse: the windows from session 1 to use for training (9:36 in
%   paper)
%   krWinToUse: the windows from session 2 to use for testing (6:54 in
%   paper)
%   preproc: the preprocessing applied to the competition data (should
%   match that used to load the krData)
%   resultRoot: directory where results should be stored
%   figTag: additional naming information for plots (e.g. '_test')
%   clusterOptions: options for cluster permutation test (see
%   clusterPermTestPooleSub_fullTime.m)
%   behav_pattern: either 'all', or a matrix describing which session 2
%   behavioral patterns to plot, e.g. [0, 1, 1, 1] or [0, 0, 1, 1;0, 1, 1,
%   1]
%   doSubZ: If 1, zscore within-subject
%   isStudy: If 1, use study trials instead of recall trials

numSub = length(subjects);

if strcmp(behav_pattern, 'all')
    behav_str = behav_pattern;
    behav_pattern = [0, 0, 0, 0; ...
        0, 0, 0, 1; ...
        0, 0, 1, 0; ...
        0, 0, 1, 1; ...
        0, 1, 0, 0; ...
        0, 1, 0, 1; ...
        0, 1, 1, 0; ...
        0, 1, 1, 1; ...
        1, 0, 0, 0; ...
        1, 0, 0, 1; ...
        1, 0, 1, 0; ...
        1, 0, 1, 1; ...
        1, 1, 0, 0; ...
        1, 1, 0, 1; ...
        1, 1, 1, 0; ...
        1, 1, 1, 1];
elseif size(behav_pattern, 1) > 1
    behav_str = num2str(reshape(behav_pattern', 1, []));
else
    behav_str = num2str(behav_pattern);
end
behav_str(isspace(behav_str)) = [];

if doSubZ
    subZ_str = '_subZ';
else
    subZ_str = '';
end

if isStudy
    study_str = '_study';
else
    study_str = '';
end

fnameString = [resultRoot '%s/KRanalysis_TGM_' preproc study_str '.mat'];
figureString = [resultRoot 'figures/KR_TGM_FirstCorrect-R4_%s_' preproc '_' study_str behav_str subZ_str figTag '.pdf'];
clusterFileString = [resultRoot 'clusters_pVals_histograms_TGM_thresh%0.2f_' preproc '_' study_str behav_str subZ_str figTag '.mat'];
clusterBootString = [resultRoot 'clusters_pVals_bootstrap_TGM_thresh%0.2f_' preproc '_' study_str behav_str subZ_str figTag '.mat'];
%%
krTGM_R = [];
krTGM_F = [];
firstRoundCorr_R = [];
firstRoundCorr_F = [];
IndividualSubjectDataR = cell(numSub, 1);
IndividualSubjectDataF = cell(numSub, 1);
IndividualSubjectFirstCorrR = cell(numSub, 1);
IndividualSubjectFirstCorrF = cell(numSub, 1);
num_good = 0;
num_tot = 0;
subjectSampleInds_R = [];
subjectSampleInds_F = [];
for i = 1:numSub
    load(sprintf(fnameString, subjects{i}));
    
    if strcmp(behav_str, 'all')
        goodItems = true(size(responseTraj, 1), 1);
    else
        goodItems = false(size(responseTraj, 1), 1);
        for i_patt = 1:size(behav_pattern, 1)
            goodItems = goodItems | ismember(responseTraj, behav_pattern(i_patt, :), 'rows');
        end
    end
    
    goodItems = goodItems & (sum(responseTraj, 2) > 1 | sum(responseTraj(:,1:3),2) > 0);
    
    num_good = num_good + sum(goodItems);
    num_tot = num_tot + length(goodItems);
    
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
    
    IndividualSubjectDataR{i} = krTrajToAdd_R;
    IndividualSubjectDataF{i} = krTrajToAdd_F;
    
    krTGM_R = cat(1, krTGM_R, krTrajToAdd_R);
    krTGM_F = cat(1, krTGM_F, krTrajToAdd_F);
    
    subjectSampleInds_R = cat(1, subjectSampleInds_R, i*ones(size(krTrajToAdd_R, 1), 1));
    subjectSampleInds_F = cat(1, subjectSampleInds_F, i*ones(size(krTrajToAdd_F, 1), 1));
end

fprintf('%d/%d trials\n', num_good, num_tot);

[~, numRounds, numKRT, numCT] = size(krTGM_R);

% Checking that everything is gaussian
ksTestOutput = nan(2, numKRT, numCT, 2);
for t1 = 1:numKRT
    for t2 = 1:numCT
        ksTestOutput(1,t1,t2, 1) = kstest(squeeze(krTGM_R(:,1:2,t1,t2)));
        ksTestOutput(1,t1,t2, 2) = kstest(squeeze(krTGM_F(:,1:2,t1,t2)));

        ksTestOutput(2,t1,t2, 1) = kstest(squeeze(krTGM_R(:,3:4,t1,t2)));
        ksTestOutput(2,t1,t2, 2) = kstest(squeeze(krTGM_F(:,3:4,t1,t2)));
    end
end

ctWinTime = -100:20:960;
ctWinTime = ctWinTime(compWinToUse);
krWinTime = winTime(krWinToUse);
krTGM_F = krTGM_F(:,:,krWinToUse, :);
krTGM_R = krTGM_R(:,:,krWinToUse, :);
numKRT = length(krWinToUse);

%% Plot desired value
numF = size(krTGM_F, 1);
numR = size(krTGM_R, 1);
numCT_short = size(krTGM_F, 4);

diffMatR = nan(numR, numKRT, numCT_short);
diffMatF = nan(numF, numKRT, numCT_short);
for i = 1:numR
    diffMatR(i,:,:) = squeeze(krTGM_R(i,firstRoundCorr_R(i),:,:) - krTGM_R(i,4,:,:));
end
for i = 1:numF
    diffMatF(i,:,:) = squeeze(krTGM_F(i,firstRoundCorr_F(i),:,:) - krTGM_F(i,4,:,:));
end

if doSubZ
    for i = 1:numSub
        numSamp_R = sum(subjectSampleInds_R == i);
        subMatR = diffMatR(subjectSampleInds_R == i, :, :);
        subMatF = diffMatF(subjectSampleInds_F == i, :, :);
        catMat = zscore(cat(1, subMatR, subMatF), 1);
        diffMatR(subjectSampleInds_R == i, :, :) = catMat(1:numSamp_R, :, :);
        diffMatF(subjectSampleInds_F == i, :, :) = catMat((numSamp_R+1):end, :, :);
    end
end


[~, rep_corr] = ttest2(diffMatR, diffMatF, 'Tail', 'right', 'Vartype', 'unequal');
rep_corr = squeeze(rep_corr);

[~, rep_opp] = ttest2(diffMatR, diffMatF, 'Tail', 'left', 'Vartype', 'unequal');
rep_opp = squeeze(rep_opp);

rep_diff = squeeze(mean(diffMatR, 1) - mean(diffMatF, 1));

rep = zeros(size(rep_opp));
rep(rep_opp < 0.10) = -1;
rep(rep_opp < 0.05) = -2;
rep(rep_opp < 0.01) = -3;
rep(rep_opp < 0.001) = -5;

rep(rep_corr < 0.10) = 1;
rep(rep_corr < 0.05) = 2;
rep(rep_corr < 0.01) = 3;
rep(rep_corr < 0.001) = 5;
%%
f1 = figure;
imagesc(krWinTime, ctWinTime, rep_diff');
ylabel('Competition Training Time (ms)');
xlabel('Swahili Learning Test Time (ms)');
colorbar
title(sprintf('Average Difference in Competition Drop\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);
set(f1, 'Color', 'w');
export_fig(f1, sprintf(figureString, 'abs'));

f2 = figure;
imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
ylabel('Competition Training Test Time (ms)');
xlabel('Swahili Learning Test Time (ms)');
colorbar
title(sprintf('P value Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);
set(f2, 'Color', 'w');
export_fig(f2, sprintf(figureString, 'PVal'));
%%
dataToCluster = cat(1, diffMatR, diffMatF);
conditionLabels = [zeros(numR, 1); ones(numF, 1)];
clusterFile = sprintf(clusterFileString, clusterOptions.pValThresh);

if exist(clusterFile, 'file')
    load(clusterFile)
else
    [clusters, pVals, permutationClusters, ...
        permutationHist, permutationSize] = clusterPermTestPooledSub_fullTime(dataToCluster, conditionLabels, clusterOptions);
    save(clusterFile, ...
        'clusters', 'pVals', 'permutationHist', 'permutationSize', ...
        'permutationClusters');
end
%%
if isempty(clusters)
    fprintf('No significant clusters found.\n')
    return
end

[~, sigClust] = min(pVals(:,2));
sizeSigClust = size(clusters{sigClust},1);
f3 = figure;
histogram(permutationSize, 'FaceAlpha', 1);
hold on
line([sizeSigClust, sizeSigClust], [0, 300]);
legend({'Permuted Cluster Sizes', 'True Max Cluster Size'});
title('Histogram of Permuted Cluster Sizes')
xlabel('Cluster size')
ylabel('Number of permuted clusters')
set(gca, 'FontSize', 18);
set(f3, 'Color', 'w');
export_fig(f3, sprintf(figureString, 'clusterSizeHist'));
%%
f4 = figure;

imagesc(krWinTime, ctWinTime, rep', [-5, 5]);
colorbar
hold on
color_text = {'ro', 'go', 'ko', 'mo', ...
        'rx', 'gx', 'kx', 'mx', ...
        'r*', 'g*', 'k*', 'm*', ...
        'rd', 'gd', 'kd', 'md',};
nr = length(krWinTime);
for i = 1:length(clusters)
    if pVals(i,2) <= 1.05
        clust = clusters{i};
        for j = 1:length(clust)
            row = mod(clust(j), nr);
            col = floor(clust(j)/nr) + 1;
            if row == 0
                row = nr;
                col = col - 1;
            end
            scatter(krWinTime(row), ctWinTime(col), color_text{i});
        end
    end
end
ylabel('Competition Training Time (ms)');
xlabel('KR Test Time (ms)');
title(sprintf('P value Difference\nRemembered vs Forgotten'));
set(gca, 'FontSize', 18);
set(f4, 'Color', 'w');

export_fig(f4, sprintf(figureString, ['PVal-Cluster_thresh' num2str(clusterOptions.pValThresh)]));
%%
clusterBootFile = sprintf(clusterBootString, clusterOptions.pValThresh);
if exist(clusterBootFile, 'file')
    load(clusterBootFile);
else
    [trueClusterT, permClusterT, bootGrid] = bootstrapCluster_KRTGM(clusters{sigClust}, ...
        clusterOptions, ...
        IndividualSubjectDataR, ...
        IndividualSubjectDataF, ...
        IndividualSubjectFirstCorrR, ...
        IndividualSubjectFirstCorrF);
    save(clusterBootFile, ...
        'trueClusterT', 'permClusterT', 'bootGrid');
end
f5 = figure;
imagesc(krWinTime, ctWinTime, bootGrid');
colorbar
caxis([0, 1])
xlabel('Swahili Learning Test Time (ms)');
ylabel('Competition Time (ms)');
title(sprintf('Proportion of Bootstrap Draws\nCluster Participation'));
set(gca, 'FontSize', 18);
set(f5, 'Color', 'w');
export_fig(f5, sprintf(figureString, 'ClusterBootstrap'));
end