% Makes stimulus set for the retrieval phase of the competition
% experiment
% Currently makes 4 blocks with 48 pairs in each block, a mix of high/low
% competition stimuli. Can easily be changed
clear;
% Seed this with the subject ID if you want to be able to reproduce
rng('shuffle');
answer=inputdlg({'Subject'},'',1,{'A'});
subject = char(answer(1,:)); %The subject the stimulus set is intended for.
if ~exist(subject, 'dir')
    mkdir(subject);
elseif exist(['./' subject '/competition.mat'], 'file')
    s = input('There already exist stimuli for this subject. Do you want to continue? y/n  ', 's');
    if s == 'n'
        return;
    end
end

% Time to display each pair
TimePerPres = 3;
% Time to display instructions
TimePerPrompt = 4;
% Time for the beginning fixation
TimePerFix = 3;
% Inter trial fixation time
ITIpres=1;

% Instructions
prompt = sprintf('Fill in the blanks\nPress the spacebar\nwhen finished.');
% Load stimuli
pairs = [importdata('./stim/C_RetC.txt'); importdata('./stim/C_RetN.txt')];
numP = length(pairs);
%Number of pairs per block
nblocks = 4;
pairs_per_block = numP/nblocks;
% Port/trigger numbers: first half competitive, second half not competitive
% port_w = 1:length(pairs);
pairsPerCat = 8;
numCat = length(pairs)/(pairsPerCat*2);
port_w = ones(1, pairsPerCat);
for p = 2:numCat
    port_w = cat(2, port_w, p*ones(1, pairsPerCat));
end
port_w = [port_w (port_w + numCat)];
% Port offset for presentation block
% Instructions = 255
% Fixation = 0;
pportStart = numCat; 
clear stories
clear s_par
clear s_len

stories = cell(nblocks);
s_par = cell(nblocks);
s_len = cell(nblocks);


% build trials
for i = 1:nblocks
    
    count = 1;
    stories{i,count} = {prompt, '+'};
    s_par{i,count} = [255, 0];
    s_len{i,count} = [TimePerPrompt, TimePerFix];
    
    rand_p = randsample(numP, pairs_per_block);
    
    for j = 1:length(rand_p)
        count = count + 1;
        stories{i, count} = {sprintf(pairs{rand_p(j)}), '+'};
        s_par{i, count} = [port_w(rand_p(j)) + pportStart, 0];
        s_len{i, count} = [TimePerPres, ITIpres];
    end
    count = count + 1;
    stories{i, count} = {'+'};
    s_par{i, count} = 0;
    s_len{i, count} = TimePerFix;
    
    
    new_ind = ones(1,numP);
    for j = 1:length(rand_p)
        new_ind = new_ind & (1:numP ~= rand_p(j));
    end
    
    pairs = pairs(new_ind);
    port_w = port_w(new_ind);
    numP = numP - pairs_per_block;
end


%Converting to single cell array
for i = 1:nblocks
    
    numTrials = size(stories, 2);
    tmp.story = cell(1, numTrials*2-1);
    tmp.parPort = [];
    tmp.len = [];
    counter = 1;
    for j = 1:numTrials
        trial = stories{i,j};
        for k = 1:length(trial)
            tmp.story{counter} = trial{k};
            counter = counter + 1;
        end
        tmp.parPort = [tmp.parPort s_par{i, j}];
        tmp.len = [tmp.len s_len{i, j}];
    end
    
    
    experiment(i).story = tmp.story; %#ok<*SAGROW>
    experiment(i).parPort = tmp.parPort;
    experiment(i).storyLength = tmp.len;
    experiment(i).storyTime = cumsum([0,tmp.len]);
end

save(['./' subject '/competition.mat'], 'experiment');

