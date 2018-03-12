function [featData, labels, winTime, options, featOptions] = ...
    preprocPipelineFreq(sub, experi, dataRoot, varargin)
% preprocPipeline: a preprocessing pipeline for CompEEG and CompEEG__KR
% data. Runs a user-selected set of preprocessing procedures (set with the
% optional struct options) and then extracts the features for use in
% classification. Note that all intermediate files are saved.
%
% Inputs:
%   sub: the subject to be processed, e.g. 'AA'
%   experi: the experiment to be processed, either 'CompEEG' or
%   'CompEEG__KR'
%   dataRoot: the location of the data. Should contain a subdir named
%   preproc-partial and further subdirectories for each subject
%   options: (optional) sets the parameters to use when preprocessing,
%   which include:
%       isVis: set to true if using visually inspected data
%       HP: the edge of the high-pass filter to apply (nan if none)
%       LP: the edge of the low-pass filter to apply (nan if none)
%       N: the center of the notch filter to apply (nan if none)
%       doRef: set to true to rereference electrodes to the group mean
%       (recommended)
%       doBase: baseline correct using prestimulus data (recommended)
%       runICA: compute ICA weights and subtract to remove blinks
%   featOptions: (optional) options for feature extraction. See
%       extractFeatures.m for details
%
% Outputs:
%   featData: the data in the form of samples x features
%   labels: the labels corresponding to these data
%   winTime: the time corresponding to each timepoint in featData (relative
%       to stimulus onset).
%   options: the basic preprocessing options used
%   featOptions: the feature extraction options used

clear EEG;
eeglab;

% Intermediate files. This script expects the originals to be in this
% directory as well.
saveInterDir = [dataRoot '/preproc-partial/' sub '/'];
% Final files
preProcFinalRoot = [dataRoot '/preproc-final/' sub '/'];
if ~exist(preProcFinalRoot, 'dir');
    mkdir(preProcFinalRoot);
end

% Epochs
epochWin = [-0.3, 2];

if nargin > 4
    options = varargin{1};
    featOptions = varargin{2};
elseif nargin > 3
    options = varargin{1};
    featOptions = struct;
else
    options = struct;
    featOptions = struct;
end

% Visual Inspection
if ~isfield(options, 'isVis')
    options.isVis = true;
end
% High Pass Filter
if ~isfield(options, 'HP')
    options.HP = 2;
end
% Low Pass Filter
if ~isfield(options, 'LP')
    options.LP = 200;
end
% Notch Filter(s)
if ~isfield(options, 'N')
    options.N = 60;
end
% Re-reference electrodes to group mean
if ~isfield(options, 'doRef')
    options.doRef = true;
end
% Limits of band to extract
if ~isfield(options, 'bpFind')
    options.bpFind = [4, 8];
end
% Name of band (for filename ease)
if ~isfield(options, 'bpName')
    options.bpName = 'theta';
end

%% EEGLab Pipeline

% Before running this pipeline, you should use the EEGLab GUI to make any
% rejections based on visual inspection. The expected filename will be of
% the form [exp '_' sub  '_Vis.set']

fnameRoot = [experi '_' sub];
if options.isVis
    fnameRoot = [fnameRoot '_Vis'];
end

% Load data
EEG = pop_loadset('filename',[fnameRoot '.set'],'filepath', saveInterDir);
EEG = eeg_checkset( EEG );
saveFname = [EEG.setname];

% HP/LP/BP Filter
BP_fname = ['_BP' num2str(options.HP) '-' num2str(options.LP)];
HP_fname = ['_HP' num2str(options.HP)];
LP_fname = ['_LP' num2str(options.LP)];
if ~isnan(options.HP) && ~isnan(options.LP) && ...
        ~exist([saveInterDir saveFname BP_fname '.set'], 'file');
    fprintf('BP filtering\n')
    EEG = pop_eegfiltnew(EEG, options.HP, options.LP, 846, 0, [], 1);
    EEG.setname=[EEG.setname BP_fname];
    saveFname = [saveFname BP_fname];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif ~isnan(options.HP) && ~isnan(options.LP) && ...
        exist([saveInterDir saveFname BP_fname '.set'], 'file');
    saveFname = [saveFname BP_fname];
    EEG = pop_loadset('filename',[saveFname '.set'],'filepath', saveInterDir);
elseif ~isnan(options.HP) && isnan(options.LP) && ...
        ~exist([saveInterDir saveFname HP_fname '.set'], 'file');
    fprintf('HP filtering\n')
    EEG = pop_eegfiltnew(EEG, options.HP, [], 846, 0, [], 1);
    EEG.setname=[EEG.setname HP_fname];
    saveFname = [saveFname HP_fname];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif ~isnan(options.HP) && isnan(options.LP) && ...
        exist([saveInterDir saveFname HP_fname '.set'], 'file')
    saveFname = [saveFname HP_fname];
    EEG = pop_loadset('filename',[saveFname '.set'],'filepath', saveInterDir);
elseif isnan(options.HP) && ~isnan(options.LP) && ...
        ~exist([saveInterDir saveFname LP_fname '.set'], 'file')
    fprintf('LP filtering\n')
    EEG = pop_eegfiltnew(EEG, [], options.LP, 846, 0, [], 1);
    EEG.setname=[EEG.setname LP_fname];
    saveFname = [saveFname LP_fname];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif isnan(options.HP) && ~isnan(options.LP) && ...
        exist([saveInterDir saveFname LP_fname '.set'], 'file')
    saveFname = [saveFname LP_fname];
    EEG = pop_loadset('filename',[saveFname '.set'],'filepath', saveInterDir);
end

% Notch Filter
if ~isnan(options.N) && ...
        ~exist([saveInterDir saveFname '_N' num2str(options.N) '.set'], 'file');
    fprintf('Notch filtering\n');
    EEG = pop_eegfiltnew(EEG, options.N - 5, options.N + 5, 846, 1, [], 1);
    EEG.setname=[EEG.setname '_N' num2str(options.N)];
    saveFname = [saveFname '_N' num2str(options.N)];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif ~isnan(options.N)
    saveFname = [saveFname '_N' num2str(options.N)];
    EEG = pop_loadset('filename',[saveFname '.set'],'filepath', saveInterDir);
end

% Re-reference the electrodes to the group mean
if options.doRef && ~exist([saveInterDir saveFname '_Ref.set'], 'file')
    fprintf('Re-referencing\n');
    EEG = pop_reref( EEG, [],'exclude',[63 65:72] );
    EEG.setname=[EEG.setname '_Ref'];
    saveFname = [saveFname '_Ref'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
elseif options.doRef
    saveFname = [saveFname '_Ref'];
    EEG = pop_loadset('filename',[saveFname '.set'],'filepath', saveInterDir);
end

% Hilbert transform
if ~any(isnan(options.bpFind)) && ~exist([saveInterDir saveFname '_Hilbert-' options.bpName '.set'], 'file')
    hilpow = hilbertTransform(EEG, options.bpFind);
    EEG.data = squeeze(hilpow);
    EEG.setname=[EEG.setname '_Hilbert-' options.bpName];
    saveFname = [saveFname '_Hilbert-' options.bpName];
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
end


% Parse into Epochs
if ~exist([saveInterDir saveFname '_Epochs.set'], 'file')
    EEG = pop_epoch( EEG, {  }, epochWin, 'newname', [EEG.setname '_Epochs'], 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    saveFname = [saveFname '_Epochs'];
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath', saveInterDir);
else
    saveFname = [saveFname '_Epochs'];
    EEG = pop_loadset('filename',[saveFname '.set'],'filepath', saveInterDir);
end

% Parse data and apply correct labels
if strcmp(experi, 'CompEEG')
    [data, labels, time] = getDataLabels_Sess1(EEG);
else
    [data, labels, time] = getDataLabels_Sess2(sub, EEG);
end

% Checking that the data/label parsing was correct
% Roughly half the trials should belong to class 1. If not using visually
% inspected data, exactly half should belong to class 1.
fprintf('%d/%d trials belong to class 1\n', sum(labels(:,1) == 1), size(labels, 1));
% If using visually inspected data and processing CompEEG__KR, make sure
% that this plot has 4 bars at 60, 180, 300, and 420
figure;
bar(labels(:,1));
title(sub);

preProcOptions = options; %#ok<*NASGU>
save([preProcFinalRoot saveFname '.mat'], 'data', 'labels', 'time', 'preProcOptions');

% Extract Relevant Features
[featData, labels, winTime, featOptions] = extractFeatures([preProcFinalRoot saveFname '.mat'], featOptions);

save([preProcFinalRoot saveFname '_Features_Overlap_Time.mat'], 'featData', 'labels', 'featOptions', 'winTime');

end