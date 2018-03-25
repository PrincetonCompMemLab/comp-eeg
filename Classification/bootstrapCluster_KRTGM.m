function [trueClusterT, permClusterT, bootGrid] = ...
    bootstrapCluster_KRTGM(clusterToUse, ...
    options, ...
    IndividualSubjectDataR, ...
    IndividualSubjectDataF, ...
    IndividualSubjectFirstCorrR, ...
    IndividualSubjectFirstCorrF)
% bootstrapCluster_KRTGM bootstraps over subjects and recomputes
% cluster test statistic for Session 3 prediction from Session 2
% competition drop. Returns the true cluster statistic and the grid of
% points showing the proportion of draws a given timepoint was significant
%
% Inputs:
%   clusterToUse: the cluster we are interested in scoring (true values)
%   options: set of options for test and cluster forming
%       maxPermutations: the number of draws to perform
%       For other options see findClusters.m
%   IndividualSubjectDataR: Cell array containing each subject's data for
%   condition 1
%   IndividualSubjectDataF: Cell array containing each subject's data for
%   condition 2
%   IndividualFirstCorrR: Cell array containing for each subject the first
%   round in which a given item was correctly recalled in condition 1
%   IndividualFirstCorrF: Cell array containing for each subject the first
%   round in which a given item was correctly recalled in condition 2
%
% Outputs:
%   trueClusterT: the score of the cluster of interest
%   permClusterT: the scores of that cluster over bootstrap draws
%   bootGrid: the grid of all timepoints showing the fraction of bootstrap
%   draws for which a given timepoint was significant
numPerms = options.maxPermutations;
numSubjects = size(IndividualSubjectDataR, 1);

% True Cluster Size
[krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F] = collectData(IndividualSubjectDataR, IndividualSubjectDataF, ...
    IndividualSubjectFirstCorrR, IndividualSubjectFirstCorrF);

trueClusterT = scoreCluster(krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F, clusterToUse, options);

permClusterT = nan(numPerms, 1);
[~, ~, numTKR, numTC] = size(krTGM_R);
bootGrid = zeros(numTKR, numTC);
subjectList = 1:numSubjects;
for p = 1:numPerms
    tic
    subjectSample = datasample(subjectList, numSubjects);
    
    [krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F] = collectData( ...
        IndividualSubjectDataR(subjectSample), IndividualSubjectDataF(subjectSample), ...
    IndividualSubjectFirstCorrR(subjectSample), IndividualSubjectFirstCorrF(subjectSample));

    [permClusterT(p), diffMatR, diffMatF] = scoreCluster(krTGM_R, krTGM_F, ...
        firstRoundCorr_R, firstRoundCorr_F, clusterToUse, options);
    
    
    dataToCluster = cat(1, diffMatR, diffMatF);
    numR = size(diffMatR, 1);
    numF = size(diffMatF, 1);
    conditionLabels = [zeros(numR, 1); ones(numF, 1)];
    clustersBoot = findClusters(dataToCluster, conditionLabels, options);
    [~, uniClustInd] = unique(cellfun(@num2str, ...
    cellfun(@(x) reshape(x, 1, []), clustersBoot, 'UniformOutput', false), ...
    'UniformOutput', false));
    clustersBoot = clustersBoot(uniClustInd);
    
    clustersBootInds = clustersBoot;
    
    for iClust = 1:length(clustersBootInds)
        bootGrid(clustersBootInds{iClust}) = bootGrid(clustersBootInds{iClust}) + 1;
    end
    
    toc
    disp(p)
end

bootGrid = bootGrid./numPerms;
end

function [krTGM_R, krTGM_F, firstRoundCorr_R, firstRoundCorr_F] = collectData(IndividualSubjectDataR, IndividualSubjectDataF, ...
    IndividualSubjectFirstCorrR, IndividualSubjectFirstCorrF)

numSubjects = length(IndividualSubjectDataR);

krTGM_R = [];
krTGM_F = [];
firstRoundCorr_R = [];
firstRoundCorr_F = [];
for s = 1:numSubjects
    krTGM_R = cat(1, krTGM_R, IndividualSubjectDataR{s});
    krTGM_F = cat(1, krTGM_F, IndividualSubjectDataF{s});
    firstRoundCorr_R = cat(1, firstRoundCorr_R, IndividualSubjectFirstCorrR{s});
    firstRoundCorr_F = cat(1, firstRoundCorr_F, IndividualSubjectFirstCorrF{s});
end

end

function [clusterT, diffMatR,diffMatF] = scoreCluster(krTGM_R, krTGM_F, ...
    firstRoundCorr_R, firstRoundCorr_F, cluster, options)

numF = size(krTGM_F, 1);
[numR, ~, numKRT, numCT] = size(krTGM_R);

diffMatR = nan(numR, numKRT, numCT);
diffMatF = nan(numF, numKRT, numCT);
for i = 1:numR
    diffMatR(i,:,:) = squeeze(krTGM_R(i,firstRoundCorr_R(i),:,:) - krTGM_R(i,4,:,:));
end
for i = 1:numF
    diffMatF(i,:,:) = squeeze(krTGM_F(i,firstRoundCorr_F(i),:,:) - krTGM_F(i,4,:,:));
end

[~, ~, ~, stats] = ttest2(diffMatR, diffMatF, 'Tail', options.tail, 'Vartype', 'unequal');

clusterT = sum(stats.tstat(cluster));
end