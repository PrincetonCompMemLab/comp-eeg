function corrAnswers = sortKRTrajBehave(sub, resultRoot, stimRoot)
% sortKRTrajBehave: gets the correctly answered session 3 items
%
% Inputs:
%   sub: subject to process, e.g. 'AA'
%   resultRoot: directory where results are be stored
%   stimRoot: directory where stimuli are stored
%
% Outputs:
%   corrAnswers: boolean indicating which items were remembered in Session
%   3

fResPrefix = [resultRoot sub '/CompEEG__KR_'];
KRanswer = importdata([stimRoot 'KR_answer.txt']);

load([resultRoot 'answers/' sub '_answers.mat']);
numQ = length(answerList);
corrAnswers = false(numQ, 1);
for a = 1:numQ
    corrAnswers(a) = strcmpi(answerList{a}, KRanswer{a});
end


save([fResPrefix sub '_corrAnswers.mat'], 'corrAnswers');

end