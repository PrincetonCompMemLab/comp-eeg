function RTs = PresentStim(par)
KbName('UnifyKeyNames');
state.esc=KbName('escape');
state.res=KbName('space');
state.atscanner = par.atscanner;
state.pportTime = par.pportTime;

eventfilename = sprintf('%s_%s_events.txt',par.subject,par.expname);
state.eventfile = fopen(eventfilename,'a');

state.bgcolor             = par.bgcolor;
state.txcolor             = par.txcolor;
state.xPixel              = 512 - par.shift;

FlushEvents;
WaitSecs(0); % preload this mex


%% Introduction
state.window = SetUpTheScreen;
state.ifi = Screen('GetFlipInterval', state.window);

try
    % Screen 1
    Screen('TextSize',state.window,64);
    Screen('FillRect',state.window,state.bgcolor);
    %     rect = Screen(state.window,'Rect');
    %     [x, y] = RectCenterd(rect);
    
    
    DrawFormattedText(state.window, sprintf('Waiting...'), state.xPixel - length(sprintf('Waiting...'))*8, 'center', state.txcolor);
    
    %     reply = GetEchoStringVert(state.window, sprintf('Meow?'), x - length(sprintf('Meow?'))*20, y-100, state.txcolor, []);
    
    Screen(state.window,'Flip');
    
    if state.atscanner,
        %         state.parport = digitalio('parallel','LPT1');
        %         hwlines = addline(state.parport,0:7,'out');
        %         putvalue(state.parport,0);
        config_io;
        outp(49200, 0);
    end
    %% fixation
    KbWait
    % Flip again just to get an experiment start time wrt the screen update
    DrawFormattedText(state.window, sprintf('Waiting...'), state.xPixel - length(sprintf('Waiting...'))*8, 'center', state.txcolor);
    vbl = Screen(state.window,'Flip');
    
    seq=1;
    RTs(1:length(par.story))=struct('stim','','indW',0,'start','length','resp',[],'key',[],'RT',[]);
    
    state.ExpStartTime = vbl;
    cur_stim.length = par.FixTime;
    cur_stim.stim = '+';
    
    cur_stim.parPortNum = par.parPort(1);
    cur_stim.pportTime = par.pportTime;
    cur_stim.stim_start = state.ExpStartTime +0.05;
    cur_stim.stim_length = par.FixTime;
    RTs(1).stim = cur_stim.stim;
    RTs(1).indW = 1;
    
    state.FixTime             = 0;
    state.NextOnsetTime = cur_stim.stim_start;
    state.drop_first=1;
    [cur_stim, state, RT] = Fix(state, cur_stim, RTs(1));
    state.drop_first=0;
    state.FixTime             = par.FixTime;
    
    
    
    
    for run = 1:length(par.story)  %% run = number of the current sentence
        
        
        seq=seq+1;
        cur_stim.stim = par.story{run};
        cur_stim.parPortNum = par.parPort(run);
        cur_stim.pportTime = par.pportTime;
        cur_stim.stim_start = par.storyTime(run);
        cur_stim.stim_length = par.storyLength(run);
        RTs(seq).stim = cur_stim.stim;
        RTs(seq).indW = run;
        
        [cur_stim, state, RTs(seq)] = ShowImgSync(state, cur_stim, RTs(seq));
        %         disp(seq);
        %         save testing.mat RTs
    end
    
    Screen(state.window,'FillRect',state.bgcolor);
    GiveBackTheScreen;
    save(sprintf('%s_%s.mat',par.subject,par.expname),'RTs');
    
    
catch %#ok<CTCH>
    GiveBackTheScreen;
    %     WriteStructsToText(filename,RTs);
    rethrow(lasterror); %#ok<LERR>
end

end

%%
%%%% PRIVATE FUNCTIONS %%%%


function [cur_stim, state, RT] = Fix(state, cur_stim, RT)
cur_stim.parPortNum =1;
cur_stim.stim = '+';
[cur_stim, state, RT] = ShowImgSync(state, cur_stim, RT);
end


function timeElapsed(message, state, cur_stim, RT)
M = sprintf('%s\t %s \t %d',message,  ...
    RT.stim,  RT.indW);
fprintf(state.eventfile,'%s\t%f\t%.6f\n',M,GetSecs-state.ExpStartTime, cur_stim.ActualOnsetTime);
end


function [cur_stim, state, RT] = ShowImgSync(state, cur_stim, RT)

textwrite = cur_stim.stim;
startStim = GetSecs;
%DrawFormattedText(state.window,textwrite,state.xPixel- length(textwrite)*8,'center',state.txcolor);
DrawFormattedText(state.window,textwrite,'center','center',state.txcolor);
% if ~strcmp(textwrite,'+'),
%     Screen(state.window, 'FillRect' ,WhiteIndex(state.window), [0 0 600 200]);
% end

FrameStartTime = state.ExpStartTime + state.FixTime + cur_stim.stim_start;
FrameStartTime = state.NextOnsetTime;
Screen('DrawingFinished', state.window);
[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos]=Screen('Flip',state.window, FrameStartTime);
cur_stim.ActualOnsetTime = VBLTimestamp;
state.NextOnsetTime = VBLTimestamp + cur_stim.stim_length - state.ifi/2;

timeElapsed('stim', state, cur_stim, RT);
if Missed >0,
    fprintf('Missed deadline!  Start stim was %.6f and deadline was %.6f\n',startStim, FrameStartTime);
    fprintf('Missed 1 %s %.6f %.6f %.6f \n',textwrite,FrameStartTime, VBLTimestamp, FrameStartTime- VBLTimestamp );
    fprintf('Missed 1 %s %.6f %.6f %.6f \n', textwrite, VBLTimestamp, FlipTimestamp, VBLTimestamp- FlipTimestamp );
end

if state.atscanner
    %%% parport command, parPortNum pportTime
    %     putvalue(state.parport,cur_stim.parPortNum);
    outp(49200,cur_stim.parPortNum);
    WaitSecs(cur_stim.pportTime);
    %     putvalue(state.parport,0);
    outp(49200,0);
else
%     fprintf('%d\n',cur_stim.parPortNum);
end
FrameReturnTime = FrameStartTime + cur_stim.stim_length/2;
keyCodes(1:256) = 0;
if state.drop_first==1
    RT = checkForInput(startStim, state, cur_stim, RT, keyCodes);
    WaitSecs(2);
end
% keyCodes(1:256) = 0; rec=0;

w_iter = 1;
% s_time = GetSecs;
rec = 0;
while (GetSecs < (FrameReturnTime))
    if sum(double(keyCodes))==0
        [keyPressed, secs, keyCodes] = KbCheck;
    elseif rec==0
        RT = checkForInput(startStim, state, cur_stim, RT, keyCodes);
        if RT.resp==1 %if they hit spacebar, stop polling
            rec=1;
        else
            keyCodes(1:256) = 0; %reset responses to look for spacebar again
        end
    end
    WaitSecs(0.0001);    % don't pound the cpu
    w_iter = w_iter +1;
    %     disp(keyPressed)
end
RT = checkForInput(startStim, state, cur_stim, RT, keyCodes);
% if RT.stim ~= '+'
%     GiveBackTheScreen;
%     fprintf('Meow\n');
% end
end



function RT = checkForInput(startStim, state, cur_stim, RT, keyCodes)
if ~sum(double(keyCodes))==0
%     fprintf('meow\n');
    charreceived = find(keyCodes==1,1);
    if charreceived == state.res
        RT.RT = GetSecs - startStim;
        %fprintf(['\n RT:' num2str(RT.RT) '\n']);
        RT.resp =-1;
        RT.key = charreceived;
        timeElapsed(strcat('R_',num2str(charreceived)), state, cur_stim, RT);
    elseif charreceived == state.esc
        GiveBackTheScreen;
        error('experimenter terminated by pressing escape');
    end
end
end

