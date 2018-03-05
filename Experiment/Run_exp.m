% Creates dialog box for selecting which experiment to run
% At Scanner? indicates whether to send triggers to machine or not
clear;
addpath expr/
KbName('UnifyKeyNames');

prompt={'Experiment Name','Subject Name', 'At Scanner? Yes=1, No=0', 'File Number',...
    'Block Number'};
def={'KR','A','0','2','1'};
title = 'SETUP EXP';
lineNo=1;

answer=inputdlg(prompt,title,lineNo,def);
exp_name = char(answer(1,:));
par.subject = char(answer(2,:));
par.atscanner = str2num(char(answer(3,:))); %#ok<ST2NM>
par.fileNum = str2num(char(answer(4,:))); %#ok<ST2NM>
par.blockNum = str2num(char(answer(5,:))); %#ok<ST2NM>
% TR = str2num(char(answer(5,:))); %#ok<ST2NM>

% Display variables
par.pportTime       = 0.1;
par.FixTime         = 5;
par.expname        = sprintf('%s_%s_%s',exp_name, num2str(par.blockNum),num2str(par.fileNum));
par.bgcolor        = 255;
par.txcolor        = 0;
par.shift           = 10;

%%
% This sets up the file naming convention and checks if you're going to
%  overwrite a subject.
filename = sprintf('%s_%s_events_redo.txt',par.subject,par.expname);
fid = fopen(filename);
while fid ~= -1
    fclose(fid);
    disp 'this file already exists ...';
    par.fileNum = input('enter new file num');
    par.expname = sprintf('%s_%s',exp_name, num2str(par.fileNum));
    filename = sprintf('%s_%s_events.txt',par.subject,par.expname);
    fid = fopen(filename);
end

%%
%Loading stimulus files for that experiment and that subject
if strcmp(exp_name, 'Comp');
    if par.fileNum == 1
        load(['./' par.subject '/presentation.mat']);
    else
        load(['./' par.subject '/competition.mat']);
    end
else
    if par.fileNum == 1
        load(['./' par.subject '/KRpres.mat']);
    else
        load(['./' par.subject '/KRtest.mat']);
    end
end
par.story = experiment(par.blockNum).story;
if isfield(experiment(par.blockNum), 'answer')
    par.answer = experiment(par.blockNum).answer;
end
par.storyLength = experiment(par.blockNum).storyLength;
par.storyTime = experiment(par.blockNum).storyTime;
par.storyTime = par.storyTime+0.5;
par.parPort = experiment(par.blockNum).parPort;

% Outputs the reaction time and response file
if strcmp(exp_name, 'Comp') || par.fileNum == 1;
    RTs=PresentStim(par);
else
    [RTs, corr] = PresentStimQ(par);
    save(['./' par.subject '/KR_results.mat'], 'RTs', 'corr');
    
end


