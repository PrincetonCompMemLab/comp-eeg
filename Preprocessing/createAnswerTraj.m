function [responseTraj] = createAnswerTraj(sub, behaveDataRoot, stimRoot)
% createAnswerTraj: Collects the answers made by subjects during the
% CompEEG__KR experiment, and forms the answer trajectory for each item.
%
% Inputs:
%   sub: the subject to be processed, e.g. 'AA'
%   behaveDataRoot: where the behavioral output from PsychToolBox is stored
%   stimRoot: the location of the stimulus files
%
% Outputs:
%   responseTraj: the 60 x 4 boolean matrix indicating whether the subject
%   answered correctly on each block for each item.

numBlocks = 4;
numStimuli = 60;

addpath ../Experiment/

eventFile = '%s/%s_KR_%d_2_events.txt';

testFid = fopen([stimRoot 'KR_test.txt']);
testStim = textscan(testFid, '%s');
testStim = testStim{1};
fclose(testFid);

ansFid = fopen([stimRoot 'KR_answer.txt']);
ansStim = textscan(ansFid, '%s');
ansStim = ansStim{1};
fclose(ansFid);

responseTraj = nan(numStimuli, numBlocks);

if exist([behaveDataRoot sub '.zip'], 'file');
    unzip([behaveDataRoot sub '.zip'], behaveDataRoot);
    delete([behaveDataRoot sub '.zip']);
end

for iblock = 1:numBlocks
    blockFid = fopen(sprintf([behaveDataRoot eventFile], sub, sub, iblock));
    events = textscan(blockFid, '%s');
    fclose(blockFid);
    events = events{1};
    
    responseInds = find(~cellfun(@isempty,regexp(events, 'R_\w*')));
    questionInds = responseInds + 1;
    
    responses = events(responseInds);
    questions = events(questionInds);
    
    for istim = 1:numStimuli
        stimInd = ismember(questions, testStim{istim});
        if any(stimInd)
            responseToCheck = responses{stimInd}(3:end);
            
            if EditDist(responseToCheck, ansStim{istim}) < 3
                responseTraj(istim, iblock) = 1;
            else
                responseTraj(istim, iblock) = 0;
            end
            % Print to check whether edit distance threshold was good
            fprintf('%s, %s = %d\n', responseToCheck, ansStim{istim}, ...
                responseTraj(istim, iblock));
        end
    end
end
save([behaveDataRoot '/' sub '/' sub '_answerTraj.mat'], 'responseTraj')

end