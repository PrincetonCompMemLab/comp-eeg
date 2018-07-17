function rt_drop_pval = plotRT(subjects, resultRoot, ...
    rtRoot, ansRoot, stimRoot)
% plotRT: plots the trajectory of reaction time (RT) over the course of
% Session 2, for trials answered correctly in session 2. Also computes the
% predictive power of reaction time drop for Session 3 behavior and returns
% the p value. Depends on results from createTGM.m.
%
% Inputs:
%   subjects: the subjects to process, e.g. {'AA', 'BB', ...}
%   resultRoot: directory where results should be stored
%   rtRoot: directory where reaction time files are stored
%   ansRoot: directory where subject Session 2 answers are stored
%   stimRoot: directory where stimuli lists are stored
% Outputs:
%   rt_drop_pval: pvalue contrasting RT drop in Session 2 for subsequently 
%   remembered and subsequently forgotten items in Session 3

numSub = length(subjects);

figureString = [resultRoot 'figures/RT_%s.pdf'];

rt_file = [rtRoot '/%s/%s_KR_%d_2.mat'];
ans_file = [ansRoot '/%s_answers.mat'];

KRanswer = importdata([stimRoot '/KR_answer.txt']);
KRtest = importdata([stimRoot '/KR_test.txt']);

fnameString = [resultRoot '/KRanalysis_TGM_%s_full.mat'];

item_names = cat(2, KRtest, KRanswer);
num_items = size(item_names, 1);
num_rounds = 4;
only_corr_rt_traj = [];
total_corr = [];
rt_drop = [];
for i_sub = 1:numSub
    sub = subjects{i_sub};
    sub_ans_load = load(sprintf(ans_file, sub));
    sub_ans = sub_ans_load.answerList;
    
    load(sprintf(fnameString, sub));
    old_good_items = sum(responseTraj, 2) > 1 | sum(responseTraj(:,1:3),2) > 0;
    
    goodItems = false(num_items, 1);
    i_old_good = 1;
    for i = 1:num_items
        if ~any(skippedItems == i)
            goodItems(i) = old_good_items(i_old_good) == 1;
            i_old_good = i_old_good + 1;
        end
    end
    
    sub_corr = nan(num_items, 1);
    sub_resp_traj = nan(num_items, num_rounds);
    sub_rt_traj = nan(num_items, num_rounds);
    sub_rt_traj_resp_corr = nan(num_items, num_rounds);
    sub_rt_drop = nan(num_items, 1);
    j_old_item = 1;
    for j_item = 1:num_items
        if ~goodItems(j_item)
            continue
        end
        sub_corr(j_item) = strcmpi(sub_ans{j_item}, KRanswer{j_item});
        first_corr = nan;
        for k_round = 1:num_rounds
            rt_load = load(sprintf(rt_file, sub, sub, k_round));
            RTs = squeeze(struct2cell(rt_load.RTs));
            
            stimuli = RTs(1,:);
            resp = RTs(5, :);
            rt = RTs(6, :);
            
            curr_item = find(ismember(stimuli, KRtest{j_item}));
            
            curr_round_corr = strcmp(resp{curr_item}, KRanswer{j_item});
            sub_rt_traj(j_item, k_round) = rt{curr_item};
            sub_resp_traj(j_item, k_round) = strcmp(resp{curr_item}, KRanswer{j_item});
            
            if sub_resp_traj(j_item, k_round)
                sub_rt_traj_resp_corr(j_item, k_round) = rt{curr_item};
            end
            
            if curr_round_corr && isnan(first_corr)
                first_corr = k_round;
            end
        end
        if isnan(first_corr)
            goodItems(j_item) = 0;
            continue
        end
        j_old_item = j_old_item + 1;
        sub_rt_drop(j_item) = sub_rt_traj(j_item, first_corr) - sub_rt_traj(j_item, num_rounds);
    end
    
    total_corr = cat(1, total_corr, sub_corr);
    rt_drop = cat(1, rt_drop, sub_rt_drop);
    only_corr_rt_traj = cat(1, only_corr_rt_traj, sub_rt_traj_resp_corr);
end

% Compute the predictive power of RT drop from first correct round to R4

rt_drop_corr = rt_drop(total_corr ==1);
rt_drop_corr(isnan(rt_drop_corr)) = [];
rt_drop_incorr = rt_drop(total_corr ==0);
rt_drop_incorr(isnan(rt_drop_incorr)) = [];

[~, rt_drop_pval] = ttest2(rt_drop_corr, rt_drop_incorr, 'Tail', 'right', 'Vartype', 'unequal');

% Plot RT Trajectory over Session 2 (for trials that were correct in
% Session 2)
rt_corr = nanmean(only_corr_rt_traj(total_corr == 1, :), 1);
rt_corr_std = nanstd(only_corr_rt_traj(total_corr == 1, :), [], 1);
rt_incorr = nanmean(only_corr_rt_traj(total_corr == 0, :), 1);
rt_incorr_std = nanstd(only_corr_rt_traj(total_corr == 0, :), [], 1);

f1 = figure;
hold on
bar(1:4, [rt_incorr; rt_corr]')
errorbar(0.86:3.86, rt_incorr, rt_incorr_std, 'r.');
errorbar(1.14:4.14, rt_corr, rt_corr_std, 'r.');
set(gca, 'XTick', 1:4);
set(gca, 'XTickLabel',{'R1', 'R2', 'R3', 'R4'});
xlabel('Round of Testing');
ylabel('Reaction Time');
legend({'Session 3 Forgotten', 'Session 3 Remembered'})
title('Average Reaction Time for Correct Session 2 Trials');
set(gca, 'fontsize', 18);
set(gcf, 'color', 'w');
legend('Remembered', 'Forgotten');
export_fig(f1, sprintf(figureString, 'traj_avg'));

end