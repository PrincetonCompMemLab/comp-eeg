function [clusters, clusterTestStat] = findClusters(data, labels, options)
% findClusters: finds significant clusters in data according to options
%
% Inputs:
%   data: a double array in one of the following two forms:
%       numSamp x dim1
%       numSamp x dim1 x dim2
%   options: a struct with any of the following fields:
%       pValThresh - uncorrected pvalue used to threshold initial pvalues 
%           to get clusters (default 0.05)
%       minClusterSize - clusters greater than this size will be considered.
%           Note that this number should be very different depending on 
%           whether your samples are time points or sensor-time points. 
%           (default 10)
%       clusterTestStatistic - summary test statistic for a cluster. Only 
%           one supported option: 'summed_t_values' the sum of the t values
%           of the samples within the cluster (default)
%       tail - tail for ttest "right" "left" or "both"
% Outputs:
%   clusters: a cell array where each entry is an array containing
%       the (linear) indices of a cluster that meets the options criteria
%   clusterTestStat: scores for each cluster

% Are we dealing with 1D series data, or 2D series?
numDim = length(size(data));

if numDim == 2
    [numSamp, numD1Points] = size(data);
elseif numDim == 3
    [numSamp, numD1Points, numD2Points] = size(data);
else
    error('Unsupported number of dimensions');
end

% Collect differences for calculating t values (for cluster
% calculation) and p values (for thresholding)
if numDim == 2
    acrossSubjectTValues = nan(numD1Points, 1);
    acrossSubjectPValues = nan(numD1Points, 1);
else
    acrossSubjectTValues = nan(numD1Points, numD2Points);
    acrossSubjectPValues = nan(numD1Points, numD2Points);
end

% Calculate t values and p values for each time or time-sensor pair across
% subjects (uncorrected)
for t = 1:numD1Points
    if numDim == 2
        [~, acrossSubjectPValues(t), ~, stats] = ttest2(data(labels==0,t), data(labels==1,t), 'Tail', options.tail, 'Vartype', 'unequal');
        acrossSubjectTValues(t) = stats.tstat;
    else
        for se = 1:numD2Points
            [~, acrossSubjectPValues(t, se), ~, stats] = ttest2(data(labels==0, t, se), data(labels==1, t, se), 'Tail', options.tail, 'Vartype', 'unequal');
            acrossSubjectTValues(t, se) = stats.tstat;
        end
    end
end

% Find the samples that survive the p value thresholding (linear indexing)

threshold_linear = acrossSubjectPValues <= options.pValThresh;

if numDim == 2
    clusters = formCluster1D(threshold_linear);
else
    clusters = formCluster2D(threshold_linear, 0);
end

% Discard clusters that are too small
numClusters = length(clusters);
thresholdedClusters = false(numClusters, 1);
for c = 1:numClusters
    if length(clusters{c}) >= options.minClusterSize
        thresholdedClusters(c) = true;
    end
end

clusters = clusters(thresholdedClusters);
numClusters = sum(thresholdedClusters);

% Compute the cluster test statistic
clusterTestStat = nan(numClusters, 1);
for c = 1:numClusters
    % Currently only one option
    switch options.clusterTestStatistic
        case 'summed_t_values'
            clusterTestStat(c) = sum(acrossSubjectTValues(clusters{c}));
    end
end

end

function clusters = formCluster1D(threshold_linear)

inds = find(threshold_linear);

% Form Initial Clusters
clusters = {};
for i = 1:length(inds)
    curr_ind = inds(i);
    if ~isempty(clusters)
        have_added = false(length(clusters), 1);
        for j = 1:length(clusters)
            clust = clusters{j};
            for k=1:length(clust)
                clust_ind = clust(k);
                if abs(curr_ind - clust_ind) == 1
                    have_added(j) = true;
                end
            end
        end
        num_touching = sum(have_added);
        if  num_touching == 1
            clusters{have_added} = cat(1, clusters{have_added}, curr_ind);
        elseif num_touching > 1
            touching_clusters = find(have_added);
            clusters{touching_clusters(1)} = cat(1, clusters{touching_clusters(1)}, curr_ind);
            clusters_to_remove = [];
            for l = 2:num_touching
                clusters_to_remove = [clusters_to_remove touching_clusters(l)]; %#ok<*AGROW>
                clusters{touching_clusters(1)} = cat(1, clusters{touching_clusters(1)}, clusters{touching_clusters(l)});
            end
            old_clusters = clusters;
            clusters = {};
            m = 1;
            for l = 1:length(old_clusters)
                if ~any(clusters_to_remove == l)
                    clusters{m} = old_clusters{l};
                    m = m + 1;
                end
            end
        else
            clusters{j+1} = curr_ind;
        end
    else
        clusters{1} = curr_ind; %#ok<*SAGROW>
    end
end

end

function clusters = formCluster2D(threshold_linear, debug)

[nr, ~] = size(threshold_linear);

inds = find(threshold_linear);

% Form Initial Clusters
clusters = {};
for i = 1:length(inds)
    curr_ind = inds(i);
    if ~isempty(clusters)
        have_added = false(length(clusters), 1);
        for j = 1:length(clusters)
            clust = clusters{j};
            for k=1:length(clust)
                clust_ind = clust(k);
                if abs(curr_ind - clust_ind) == 1 || abs(curr_ind - clust_ind) == nr
                    if ~(mod(curr_ind, nr) == 0 && mod(clust_ind, nr) == 1) && ...
                            ~(mod(curr_ind, nr) == 1 && mod(clust_ind, nr) == 0)
                        have_added(j) = true;
                    end
                end
            end
        end
        num_touching = sum(have_added);
        if  num_touching == 1
            clusters{have_added} = cat(1, clusters{have_added}, curr_ind);
        elseif num_touching > 1
            touching_clusters = find(have_added);
            clusters{touching_clusters(1)} = cat(1, clusters{touching_clusters(1)}, curr_ind);
            clusters_to_remove = [];
            for l = 2:num_touching
                clusters_to_remove = [clusters_to_remove touching_clusters(l)]; %#ok<*AGROW>
                clusters{touching_clusters(1)} = cat(1, clusters{touching_clusters(1)}, clusters{touching_clusters(l)});
            end
            old_clusters = clusters;
            clusters = {};
            m = 1;
            for l = 1:length(old_clusters)
                if ~any(clusters_to_remove == l)
                    clusters{m} = old_clusters{l};
                    m = m + 1;
                end
            end
        else
            clusters{j+1} = curr_ind;
        end
    else
        clusters{1} = curr_ind; %#ok<*SAGROW>
    end
end

if debug
    imagesc(threshold_linear)
    hold on
    color_text = {'ro', 'go', 'bo', 'ko', 'mo', ...
        'rx', 'gx', 'bx', 'kx', 'mx', ...
        'r*', 'g*', 'b*', 'k*', 'm*', ...
        'rd', 'gd', 'bd', 'kd', 'md',};
    for i = 1:length(clusters)
        clust = clusters{i};
        for j = 1:length(clust)
            row = mod(clust(j), nr);
            col = floor(clust(j)/nr) + 1;
            if row == 0
                row = nr;
                col = col - 1;
            end
            scatter(col, row, color_text{mod(i, length(color_text)) + 1});
        end
    end
end
end