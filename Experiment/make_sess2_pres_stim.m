% Makes stimulus set for the study/presentation phase of the Swahili exp
clear;
% Seed this with the subject ID if you want to be able to reproduce
rng('shuffle');
answer=inputdlg({'Subject'},'',1,{'A'});
subject = char(answer(1,:)); %The subject the stimulus set is intended for.
if ~exist(subject, 'dir')
    mkdir(subject);
elseif exist(['./' subject '/KRpres.mat'], 'file')
    s = input('There already exist stimuli for this subject. Do you want to continue? y/n  ', 's');
    if s == 'n'
        return;
    end
end

% Time to display each pair
TimePerPres = 2;
% Time to display instructions
TimePerPrompt = 4;
% Time for the beginning fixation
TimePerFix = 3;
% Inter trial fixation time
ITIpres=1;

% Instructions
prompt = sprintf('Try to memorize the\n following word pairs.');
% Study stimuli (same order as KR_test.txt)
pairs = importdata('./stim/KR_study_ordered.txt');
% Number of stimuli
numP = length(pairs);
% Port numbers
port_w = 1:length(pairs);

% Number of blocks
nblocks= 4;
% Port/trigger number offset for instructions
% Instructions = 255
% Fixation = 0;
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
    
    rand_w = randperm(numP);
    for j = 1:length(rand_w)
        count = count + 1;
        stories{i, count} = {sprintf(pairs{rand_w(j)}), '+'};
        s_par{i, count} = [port_w(rand_w(j)), 0];
        s_len{i, count} = [TimePerPres, ITIpres];
    end
    count = count + 1;
    stories{i, count} = {'+'};
    s_par{i, count} = 0;
    s_len{i, count} = TimePerFix;
    
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

save(['./' subject '/KRpres.mat'], 'experiment');

