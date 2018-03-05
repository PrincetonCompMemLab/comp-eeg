% Makes stimulus set for the testing/retrieval phase of Swahili exp
clear;
% Seed this with the subject ID if you want to be able to reproduce
rng('shuffle');
answer=inputdlg({'Subject'},'',1,{'A'});
subject = char(answer(1,:)); %The subject the stimulus set is intended for.
if ~exist(subject, 'dir')
    mkdir(subject);
elseif exist(['./' subject '/KRtest.mat'], 'file')
    s = input('There already exist stimuli for this subject. Do you want to continue? y/n  ', 's');
    if s == 'n'
        return;
    end
end

% Minimum time to present the stimulus
TimePerPres = 2;
% Instruction presentation time
TimePerPrompt = 4;
% Beginning fixation time
TimePerFix = 3;
% Inter trial fixation interval length
ITIpres=1;

% Instructions
prompt = sprintf('Type the English\ntranslation,\nthen hit enter.');
% Swahili words
swa = importdata('./stim/KR_test.txt');
% English words
eng = importdata('./stim/KR_answer.txt');
% Number of stimuli
numP = length(swa);
% Port/trigger number assignments (determined by word order in KR_test.txt)
port_w = 1:length(swa);

% Number of blocks
nblocks= 4;
% Instructions = 255
% Fixation = 0
clear stories
clear s_par
clear s_len

% Contains the prompt stimuli (swahili words), to be shown
stories = cell(nblocks);
% Contains answer stimuli (english words), not shown to subject
answers = cell(nblocks);
% Port numbers
s_par = cell(nblocks);
% Time to display each stimulus
s_len = cell(nblocks);


% build trials
for i = 1:nblocks
    
    count = 1;
    stories{i,count} = {prompt, '+'};
    answers{i,count} = {'null', 'null'};
    s_par{i,count} = [255, 0];
    s_len{i,count} = [TimePerPrompt, TimePerFix];
    
    rand_w = randperm(numP);
    for j = 1:length(rand_w)
        count = count + 1;
        stories{i, count} = {sprintf(swa{rand_w(j)}), '+'};
        answers{i, count} = {eng{rand_w(j)}, 'null'};
        s_par{i, count} = [port_w(rand_w(j)), 0];
        s_len{i, count} = [TimePerPres, ITIpres];
    end
    count = count + 1;
    stories{i, count} = {'+'};
    answers{i, count} = {'null'};
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
        answer = answers{i,j};
        for k = 1:length(trial)
            tmp.story{counter} = trial{k};
            tmp.ans{counter} = answer{k};
            counter = counter + 1;
        end
        tmp.parPort = [tmp.parPort s_par{i, j}];
        tmp.len = [tmp.len s_len{i, j}];
    end
    
    
    experiment(i).story = tmp.story; %#ok<*SAGROW>
    experiment(i).answer = tmp.ans;
    experiment(i).parPort = tmp.parPort;
    experiment(i).storyLength = tmp.len;
    experiment(i).storyTime = cumsum([0,tmp.len]);
end

save(['./' subject '/KRtest.mat'], 'experiment');

