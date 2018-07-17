function plotBehaveResults(subjects, resultRoot, ...
    behaveDataRoot, stimRoot, figTag)
% plotBehaveResults plots summary statistics of the behavioral data from
% Session 2 and Session 3
%
% Inputs:
%   subjects: subjects to process (e.g. {'AA', 'BB', ...})
%   resultRoot: dir to save/load results (assumes subject subdirs)
%   behaveDataRoot: dir containing subject behavioral data (assumes subject 
%   subdirs)
%   stimRoot: dir where stimulus data is stored
%   figTag: additional naming information for plots (e.g. '_test')
numSub = length(subjects);

meanRespTraj = [];
meanRespTrajCorr = [];
meanRespTrajInc = [];
meanCorr = [];
krTrajList = cell(numSub, 1);
krLabelList = cell(numSub, 1);
perSubRespTrajCorr = nan(numSub, 4);
perSubRespTrajInc = nan(numSub, 4);
perSubRespTraj = nan(numSub, 4);
perSubCorr = nan(numSub, 1);
configurationPerf = zeros(16,2);
for s = 1:numSub
    load(sprintf('%s%s/%s_answerTraj.mat', behaveDataRoot, subjects{s}, subjects{s}));
    corrAnswers = sortKRTrajBehave(subjects{s}, resultRoot, stimRoot);
    
    numTraj = size(responseTraj, 1);
    noNaN = false(numTraj, 1);
    for t = 1:numTraj
        if ~any(isnan(responseTraj(t,:)))
            noNaN(t) = true;
        end
    end
    
    krTrajList{s} = responseTraj(noNaN, :);
    krLabelList{s} = corrAnswers(noNaN);
    
    for it = 1:size(krTrajList{s},1)
        if krTrajList{s}(it,1) == 0
            if krTrajList{s}(it,2) == 0
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(1,krLabelList{s}(it,1)+1) = configurationPerf(1,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(2,krLabelList{s}(it,1)+1) = configurationPerf(2,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(3,krLabelList{s}(it,1)+1) = configurationPerf(3,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(4,krLabelList{s}(it,1)+1) = configurationPerf(4,krLabelList{s}(it,1)+1)+1;
                    end
                end
            else
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(5,krLabelList{s}(it,1)+1) = configurationPerf(5,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(6,krLabelList{s}(it,1)+1) = configurationPerf(6,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(7,krLabelList{s}(it,1)+1) = configurationPerf(7,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(8,krLabelList{s}(it,1)+1) = configurationPerf(8,krLabelList{s}(it,1)+1)+1;
                    end
                end
            end
        else
            if krTrajList{s}(it,2) == 0
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(9,krLabelList{s}(it,1)+1) = configurationPerf(9,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(10,krLabelList{s}(it,1)+1) = configurationPerf(10,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(11,krLabelList{s}(it,1)+1) = configurationPerf(11,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(12,krLabelList{s}(it,1)+1) = configurationPerf(12,krLabelList{s}(it,1)+1)+1;
                    end
                end
            else
                if krTrajList{s}(it,3) == 0
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(13,krLabelList{s}(it,1)+1) = configurationPerf(13,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(14,krLabelList{s}(it,1)+1) = configurationPerf(14,krLabelList{s}(it,1)+1)+1;
                    end
                else
                    if krTrajList{s}(it,4) == 0
                        configurationPerf(15,krLabelList{s}(it,1)+1) = configurationPerf(15,krLabelList{s}(it,1)+1)+1;
                    else
                        configurationPerf(16,krLabelList{s}(it,1)+1) = configurationPerf(16,krLabelList{s}(it,1)+1)+1;
                    end
                end
            end
        end
    end
    
    meanRespTraj = cat(1, meanRespTraj, responseTraj);
    meanRespTrajCorr = cat(1, meanRespTrajCorr, responseTraj(corrAnswers,:));
    meanRespTrajInc = cat(1, meanRespTrajInc, responseTraj(~corrAnswers,:));
    meanCorr = cat(1, meanCorr, corrAnswers);
    perSubRespTrajCorr(s, :) = nanmean(responseTraj(corrAnswers, :), 1);
    perSubRespTrajInc(s, :) = nanmean(responseTraj(~corrAnswers, :), 1);
    perSubRespTraj(s, :) = nanmean(responseTraj, 1);
    perSubCorr(s) = nanmean(corrAnswers);
end
%% Configuration plots

% Total prevalence of each configuration
f = figure;
bar(sum(configurationPerf, 2))
xlim([0.5, 16.5])
set(gca, 'xtick', 1:16);
set(gca, 'xticklabels', {'0000', '0001', '0010', '0011', '0100', '0101', ...
    '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111'});
title('Prevalence of Response Configurations During Session 2');
xlabel('Response Configuration from Session 2')
ylabel('Number of Items')
set(gca, 'FontSize', 16);
set(f, 'Color', 'w');
set (f, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f, [resultRoot 'figures/Behavioral/configPrev' figTag '.pdf']);
f = figure;
bar(configurationPerf)
xlim([0, 17])
set(gca, 'xtick', 1:16);
set(gca, 'xticklabels', {'0000', '0001', '0010', '0011', '0100', '0101', ...
    '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111'});
legend({'Session 3 Forgotten', 'Session 3 Remembered'});
title(sprintf('Prevalence of Response Configurations During Session 2\nSplit by Session 3 Performance'));
xlabel('Response Configuration from Session 2')
ylabel('Number of Items')
set(gca, 'FontSize', 20);
set(f, 'Color', 'w');
set (f, 'Units', 'normalized', 'Position', [0,0,1,1]);
export_fig(f, [resultRoot 'figures/Behavioral/configPrevSplit' figTag '.pdf']);
%%
f1 = figure;
totalMean = [nanmean(perSubRespTraj, 1) nanmean(perSubCorr)];
totalStd = [nanstd(perSubRespTraj) nanstd(perSubCorr)]./sqrt(numSub);
bar(1:5, totalMean);
hold on
set(gca, 'XTickLabel',{'R1', 'R2', 'R3', 'R4', 'T'});
errorbar(1:5, totalMean, totalStd, 'r.');
xlabel('Round of Testing');
ylabel('Fraction Correct');
title(sprintf('Average Performance Across Subjects\nin Each Session 2 Round and in Session 3'));
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f1, [resultRoot 'figures/Behavioral/avgPerf_pooled' figTag '.pdf']);

f15 = figure;
h = histogram(perSubCorr, 'BinEdges', 0:0.1:1);
set(h, 'FaceAlpha', 1);
xlabel('Average Performance in Session 3')
ylabel('Number of Subjects')
title(sprintf('Histogram of Individual Subject\nPerformance in Session 3'))
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f15, [resultRoot 'figures/Behavioral/subHist_Session3' figTag '.pdf']);

f2 = figure;
corrIncMean = [nanmean(perSubRespTrajInc, 1); nanmean(perSubRespTrajCorr, 1)];
bar(1:4, corrIncMean');
hold on
set(gca, 'XTickLabel',{'R1', 'R2', 'R3', 'R4'});
errorbar(1.145:4.145, corrIncMean(2,:), nanstd(perSubRespTrajCorr)./sqrt(numSub), 'k.');
errorbar(0.855:3.855, corrIncMean(1,:), nanstd(perSubRespTrajInc)./sqrt(numSub), 'r.');
legend('Session 3 Forgotten', 'Session 3 Remembered', 'Location', 'NorthWest');
xlabel('Round of Testing');
ylabel('Fraction Correct in Session 2');
title(sprintf('Average Session 2 Performance Across Subjects\nGrouped by Session 3 Performance'));
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
export_fig(f2, [resultRoot 'figures/Behavioral/avgPerf_div_pooled' figTag '.pdf']);
end