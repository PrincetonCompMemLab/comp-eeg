function [data, labels, time] = getDataLabels_Sess1(EEG)
% getDataLabels_Sess1: converts the eeglab data structure EEG into a more
% workable format. Designed for CompEEG session 1.
%
% Inputs:
%   EEG: the eeglab data structure to be converted
%
% Outputs:
%   data: the data as a double array of size samples x channels x time
%   labels: the labels for the data
%   time: a vector corresponding the time dimension of data with
%   stimulus-relative timepoints

% Extract the epochs
epochInfo = EEG.epoch;
numEpochs = length(epochInfo);

labelInds = false(numEpochs, 1);
stimLabels = nan(numEpochs, 2); %Comp/nonComp binary, category identity

% Trigger start number known from the experiment-running code
labelStart = 19;

for i = 1:numEpochs
    if iscell(epochInfo(i).eventtype)
        epochInfo(i).eventtype = epochInfo(i).eventtype{1};
    end
    if isstr(epochInfo(i).eventtype) %#ok<*DISSTR>
        label = str2num(epochInfo(i).eventtype);%#ok<*ST2NM>
    else
        label = epochInfo(i).eventtype;
    end
    if (label > labelStart) && (label ~= 255)
        labelInds(i) = true;
        stimLabels(i, 2) = label;
        
        if label > labelStart*2
            stimLabels(i, 1) = 0;
        else
            stimLabels(i, 1) = 1;
        end
    end
end

labels = stimLabels(labelInds, :);
data = EEG.data(:,:,labelInds);
data = permute(data, [3 1 2]);
numChan = size(data, 2);
data = data(:, 1:(numChan - 8), :);
time = EEG.times;

end
