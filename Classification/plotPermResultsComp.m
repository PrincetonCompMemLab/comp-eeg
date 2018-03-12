function plotPermResultsComp(threshold, ...
    preproc, numPerm, numFolds, doZ, resultRoot)
% plotPermResultsComp: plots results of competition classification and
% permutation test - corrects for multiple comparisons and combines data
% across subjects. 
%
% Inputs:
%   threshold: the p value threshold (will be Bonferroni corrected)
%   preproc: string indicating what preprocessing has been applied to the
%   data (theta or voltage)
%   numPerm: the number of permutation draws
%   numFolds: the number of folds
%   doZ: whether to zscore or not (1 in paper)
%   resultRoot: where you want the result saved. Expecting subdir for each 
%   subject
addpath /Users/nrafidi/Documents/MATLAB/Toolboxes/export_fig/


fileString = [resultRoot 'CompEEG_' num2str(numFolds) ...
    'FCV_win*_permAccs' num2str(numPerm) ...
    '_Z' num2str(doZ) '_' preproc '.mat'];

filesToUse = dir(fileString);

%%
numFiles = length(filesToUse);

trueAccs = nan(numFiles, 1);
trueSubAccsPooled = [];
pVals = nan(numFiles, 1);
timeToPlot = nan(numFiles, 1);
for f = 1:numFiles
    load([resultRoot filesToUse(f).name]);
    try
        timeToPlot(f) = winTime(str2double(filesToUse(f).name(17:18)));
    catch
        timeToPlot(f) = winTime(str2double(filesToUse(f).name(17)));
    end
    trueAccs(f) = mean(mean(trueSubAccs));
    trueSubAccs = mean(trueSubAccs, 2);
    trueSubAccsPooled = cat(2, trueSubAccsPooled, trueSubAccs);
    permAccs = squeeze(mean(permSubAccs, 3));
    sum_perms_greater = sum(permAccs > repmat(trueSubAccs, 1, numPerm), 2);
    subPvals = (sum_perms_greater + 1)/numPerm;
    subProbs = norminv(subPvals-eps, 0, 1);
    [~, pVals(f)] = ttest(subProbs);
end

erpSubAvg = nanmean(trueSubAccsPooled);
erpSubStd = nanstd(trueSubAccsPooled);

[timeToPlot, sortInds] = sort(timeToPlot);

erpSubAvg = erpSubAvg(sortInds);
erpSubStd = erpSubStd(sortInds);
pVals = pVals(sortInds);
trueAccs = trueAccs(sortInds);

subUpp = erpSubAvg + erpSubStd;
subLow = erpSubAvg - erpSubStd;
X = [timeToPlot', fliplr(timeToPlot')];
Y = [subLow, fliplr(subUpp)];

compInds = find(subLow > 0.5);
disp(compInds)

h = figure;
hold on
fh = fill(X, Y, 'c');
set(fh, 'EdgeAlpha', 0);
plot(timeToPlot, erpSubAvg, 'b');
threshold = threshold/length(timeToPlot);
disp(threshold)
for f = 1:(numFiles-1)
    if pVals(f) <= threshold
        scatter(timeToPlot(f), trueAccs(f)+0.005, 'r*');
    end
end
legend({'Standard Deviation', 'Accuracy', 'Bonferroni Significant'});
line([min(timeToPlot), max(timeToPlot)], [0.5, 0.5], 'LineStyle', '--', 'Color', 'k');
xlim([min(timeToPlot), max(timeToPlot)]);
ylim([0.45, 0.65]);
xlabel('Time relative to stimulus onset (ms)');
ylabel('Accuracy');

title(sprintf('Competition Decoding Accuracy\nOver Time in Session 1\n%s', preproc));
set(gca, 'FontSize', 18)
set(h, 'Color', 'w');
export_fig(h, [resultRoot 'figures/compPerm_' num2str(numFolds) ...
    'FCV_permAccs' num2str(numPerm) ...
    '_Z' num2str(doZ) '_' preproc '_thresh' num2str(threshold) '.pdf']);