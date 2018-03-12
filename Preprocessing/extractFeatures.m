function [featData, labels, erpWin, options] = extractFeatures(loadFname, varargin) %#ok<*STOUT>
% extractFeatures: a function for extracting features from EEG data
% recorded in the file loadFname, with optional input parameters. Note that
% this is intended to be called from preprocPipeline.m
%
% Inputs:
%   loadFname: the filename containing the preprocessed EEG data. It should
%   contain a variable data that is samples x channels x time, a variable
%   labels, and a variable time
%   options: (optional) sets the parameters to use when extracting, which
%   include:
%       erpWinSize: the size of window to average
%       overLap: overlapping erp windows if true
%       minTime: the minimum time (relative to onset) to consider for
%       extraction.
%       maxTime: the maximum amount of time post-onset to consider for
%       extraction.
%
% Outputs:
%   featData: the data in the form of samples x features
%   labels: the labels corresponding to these data
%   erpWin: the start times of the erp windows
%   options: the feature extraction options used

% Parameters
if nargin > 1
    options = varargin{1};
else
    options = struct;
end
% Window size to average voltages (ms)
if ~isfield(options, 'erpWinSize')
    options.erpWinSize = 50;
end
% Overlap windows?
if ~isfield(options, 'overLap')
    options.overLap = true;
end
% Maximum time post-stim to consider (ms)
if ~isfield(options, 'maxTime')
    options.maxTime = 1000;
end

load(loadFname);

% Minimum time relative to stim to consider (ms)
if ~isfield(options, 'minTime')
    minTime = min(time(time >= -100));
else
    minTime = min(time(time >= options.minTime));
end

if options.overLap
    erpWin = minTime:20:options.maxTime;
else
    erpWin = minTime:options.erpWinSize:options.maxTime;
end

numErp = length(erpWin);

featData = [];

for e = 2:numErp
    erp = (time >= erpWin(e-1)) & (time < (erpWin(e-1)+options.erpWinSize));
    newData = squeeze(mean(data(:,:,erp), 3));
    featData = cat(2, featData, newData);
end

end