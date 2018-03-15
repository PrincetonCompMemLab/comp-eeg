function [clusters, monteCarloPvals, permutationClusters, ...
    permutationHist, sizePermClusters] = ...
    clusterPermTestPooledSub_fullTime(pooledData, varargin)
% clusterPermTestPooledSub_fullTime: runs the cluster permutation test 
% described in Maris & Oostenveld 2007 for two conditions in a 
% within-subjects E/MEG study
%
% Inputs:
% pooledData: a double array in one of the following two forms:
%       numSamp x dim1 x condition
%       numSamp x dim1 x dim2 x condition
%
% dim1 and dim2 could be time, space, whatever
%
% Optional Inputs:
% options: a struct with any of the following fields:
%   pValThresh - uncorrected pvalue used to threshold initial pvalues to get
%       clusters (default 0.05)
%   minClusterSize - clusters greater than this size will be considered.
%       Note that this number should be very different depending on whether
%       your samples are time points or sensor-time points. (default 10)
%   clusterTestStatistic - summary test statistic for a cluster. Only one
%       supported option: 'summed_t_values' the sum of the t values of the
%       samples within the cluster (default)
%   maxPermutations - maximum number of condition assignment permutations 
%       to try (default 1000)
%
% Outputs:
%   clusters: a cell array where each entry is an array containing
%       the (linear) indices of a cluster that meets the options criteria
%   monteCarloPvals: the Monte Carlo p values for each of the significant
%       clusters.
%   permutationClusters: the set of clusters recovered in each permutation
%   permutationHist: the set of cluster scores recovered over permutations
%   sizePermClusters: permutationHist but containing the size of the
%       clusters instead of the score specified by options

% Checking Arguments
if nargin < 2
    options = struct();
elseif nargin == 2
    options = varargin{1};
else
    error('Too many input arguments');
end

% Are we dealing with 1D series data, or 2D series?
numDim = length(size(pooledData));

if numDim == 3
    [numSamp, numD1Points, numCond] = size(pooledData);
elseif numDim == 4
    [numSamp, numD1Points, numD2Points, numCond] = size(pooledData);
else
    error('Unsupported number of dimensions');
end

% Checking the number of conditions (right now can only handle two)
if numCond ~= 2
    error('Unsupported number of conditions');
end

% Setting options to defaults if they are not provided
if ~isfield(options, 'pValThresh')
    options.pValThresh = 0.05;
end
if ~isfield(options, 'minClusterSize')
    options.minClusterSize = 10;
end
if ~isfield(options, 'clusterTestStatistic')
    options.clusterTestStatistic = 'summed_t_values';
end
if ~isfield(options, 'maxPermutations')
    options.maxPermutations = 1000;
end

[clusters, clusterTestStat] = findClusters(pooledData, options);
numClusters = length(clusters);

% Get the permutations of condition assignments
numPermutations = options.maxPermutations;
permutationClusters = cell(numPermutations, 1);
permutationTestStat = [];
numPermClusters = [];
sizePermClusters = [];
permutationHist = [];
percFlipped = nan(numPermutations, 1);

% Compute the cluster statistic on the max cluster data for each of the
% permuted condition assignments
for p = 1:numPermutations
    disp(p)
    permuted = randi([0,1], 1, numSamp);
    percFlipped(p) = 1 - nnz(permuted)/numSamp;
    permuted = [permuted + 1; 2 - permuted];
    
    permutedData = nan(size(pooledData));
    
    for s = 1:numSamp
        if numDim == 3
            permutedData(s,:,:) = pooledData(s,:,permuted(:,s));
        elseif numDim == 4
            permutedData(s,:,:,:) = pooledData(s,:,:,permuted(:,s));
        end
    end
    
    [permutationClusters{p}, permutationTestStat_single] = ...
        findClusters(permutedData, options);
    
    if isempty(permutationClusters{p})
        permutationTestStat_single = 0;
        sizePermClusters = cat(1, sizePermClusters, 0);
    else
        sizePermClusters = cat(1, sizePermClusters, max(cellfun(@(x) size(x, 1), permutationClusters{p})));
    end
    
    [maxAbs, indMax] = max(abs(permutationTestStat_single));
    permutationHist = cat(1, permutationHist, permutationTestStat_single(indMax));
    permutationTestStat = cat(1, permutationTestStat, maxAbs);
    numPermClusters = cat(1, numPermClusters, length(permutationClusters{p}));
end

% For each cluster, its MC p value is the fraction of permutation values
% that are greater than its own. This is a two sided test, so we compare
% absolute values
monteCarloPvals = nan(numClusters, 2);
for c = 1:numClusters
    monteCarloPvals(c,1) = sum(abs(permutationTestStat) > abs(clusterTestStat(c)))/length(permutationTestStat);
    monteCarloPvals(c,2) = sum(sizePermClusters > size(clusters{c},1))/length(sizePermClusters);
end
disp(nnz(numPermClusters))
disp(sum(numPermClusters > numClusters));

end