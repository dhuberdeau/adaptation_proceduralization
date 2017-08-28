function [trial_triplet, stats] = nearest_neighbor_analysis_fnc(category_data, reach_data, trials_data, probe_data, trials_a, trials_b, sub_exl)
% trial_triplet = nearest_neighbor_analysis_fnc(category_data, reach_data, trials_data, probe_data, trials_a, trials_b, sub_exl)
% 
% dipslay probe trials and their nearest neighbors for first half and
% second half of adaptation windows.
% 
% INPUTS:  N = number of trials, S = number of subjects
%   category_data - NxS matrix of the category for each trial
%   pt_data - NxS matrix of the PT latency of each trial
%   reach_data - NxS matrix of the reach direction computed for each trial
%   trials_data - NxS matrix of the trial numbers, trivial, all S same
%   probe_data - NxS matrix of boolean values whether trial is cat6 (bp)
%   
%   
% 
% OUTPUTS: cell array, trial_triplet, containing triplets of trials for each
% probe category type 3xRxS, R = repeats of the triplet in a block.
% cell array, stats, contains p-values and t-test t-values for cat2 vs.
% cat4 (rate and asymptote), and cat2 vs. cat6 (rate and asymptote)...
% {p4_rate, p4_asym, p6_rate, p6_asym; stat24a, stat24b, stat26a, stat26b};

% load('data_E1_ptFix_2TH_2cm_vManRev.mat')

% subject_exclusions = [2 5 10 11 16]; %identified as bad apples08/12/2016
subject_exclusions = sub_exl;

all_subs = 1:size(category_data,2);
keep_subs = setdiff(all_subs, subject_exclusions);
category_temp = category_data;
reach_temp = reach_data;
trials_temp = trials_data;
probe_temp = probe_data;
clear category_data pt_data reach_data trials_data

category_data = category_temp(:, keep_subs);
reach_data = reach_temp(:, keep_subs);
trials_data = trials_temp(:, keep_subs);
probe_data = probe_temp(:, keep_subs);
probe_off = ([zeros(1, size(probe_data,2)); diff(probe_data)] < 0);
clear *_temp

%% Find nearest neighbors:
%rotation onset trial 280
%rotation offset trial 400

additional_subjs_exclude = [];

n_cat4 = nan(1, size(category_data,2));
for i_sub = 1:size(category_data,2)
    n_cat4(i_sub) = sum(category_data(trials_data(:, i_sub) > trials_a(1) & trials_data(:, i_sub) <= trials_a(2), i_sub) == 4);
end
n_cat6 = nan(1, size(category_data,2));
for i_sub = 1:size(category_data,2)
    cond1 = trials_data(:, i_sub) > trials_a(1) & trials_data(:, i_sub) <= trials_a(2);
    n_cat6(i_sub) = sum(category_data(cond1, i_sub) == 2 & probe_data(cond1, i_sub));
end
n_cat4_mode = mode(n_cat4);
n_cat6_mode = mode(n_cat6);

trial_triplet_cat4_a = nan(3, n_cat4_mode, size(category_data, 2));
trial_triplet_cat6_a = nan(3, n_cat6_mode, size(category_data, 2));

for i_sub = 1:size(category_data,2)
    subset1 = trials_data(:, i_sub) > trials_a(1) & trials_data(:, i_sub) <= trials_a(2);
    reach_subset = reach_data(subset1, i_sub);
    reach_subset(reach_subset > 50 | reach_subset < -50) = nan;
    cat_subset = category_data(subset1, i_sub);
    probe_subset = probe_data(subset1, i_sub);
    probe_off_subset = probe_off(subset1, i_sub);
    cat2_inds = find(cat_subset == 2 & probe_subset == 0);
    cat4_inds = find(cat_subset == 4);
    cat6_inds = find(cat_subset == 2 & probe_subset == 1);
    cat2_p_inds = find(cat_subset == 2 & probe_off_subset == 1); % indicies of LPT trials that directly follow Aftereffect trials
    for i_4 = 1:size(trial_triplet_cat4_a, 2)
        trial_triplet_cat4_a(2, i_4, i_sub) = reach_subset(cat4_inds(i_4));
        k_before = max(cat2_inds(cat2_inds < cat4_inds(i_4)));
        k_after = min(cat2_inds(cat2_inds > cat4_inds(i_4)));
        if ~isempty(cat6_inds)
            if ismember(k_before, cat2_p_inds)
                % k_before is a post-aftereffect trial, omit it
                k_before = nan;
            end
            if ismember(k_after, cat2_p_inds)
                % k_after is a post-aftereffect trial, omit it
                k_after = nan;
            end
        end
        try
            if ~isnan(k_before)
                trial_triplet_cat4_a(1, i_4, i_sub) = reach_subset(k_before);
            end
        catch
        end
        try
            if ~isnan(k_after)
                trial_triplet_cat4_a(3, i_4, i_sub) = reach_subset(k_after);
            end
        catch boo
            warning(boo.message)
        end
    end
    if ~isempty(cat6_inds)
        for i_6 = 1:size(trial_triplet_cat6_a, 2)
            trial_triplet_cat6_a(2, i_6, i_sub) = reach_subset(cat6_inds(i_6));
            k_before = max(cat2_inds(cat2_inds < cat6_inds(i_6)));
            cat2_inds_after = cat2_inds(cat2_inds > cat6_inds(i_6));
            [k_after, inds_k_after] = min(cat2_inds_after);
            trial_triplet_cat6_a(1, i_6, i_sub) = reach_subset(k_before);
            try
                trial_triplet_cat6_a(3, i_6, i_sub) = reach_subset(cat2_inds_after(2));
            catch

            end
        end
    end
end

% trial_triplet_cat4_a = trial_triplet_cat4_a(:, 2:end, :);
figure; subplot(2,2,1); hold on;
% trial_triplet_cat4_a = trial_triplet_cat4_a./repmat(trial_triplet_cat4_a(1,:,:), 3,1,1);
%exclude subject 9: bc they are outlier
for i_sub_x = 1:length(additional_subjs_exclude)
    trial_triplet_cat4_a(:,:,additional_subjs_exclude(i_sub_x)) = nan;
end
for i_sub = 1:size(trial_triplet_cat4_a,3)
    compl_data = sum(sum(~isnan(trial_triplet_cat4_a(:,:,i_sub)),2)>0);
    if compl_data >= 3
        plot(1:3, nanmean(trial_triplet_cat4_a(:,:, i_sub), 2), '.-', 'Color', [.5 .5 .5]);
    end
end
trial_triplet_cat4_a = reshape(trial_triplet_cat4_a, 3, size(trial_triplet_cat4_a,2)*size(trial_triplet_cat4_a,3));
% trial_triplet_cat4_a = trial_triplet_cat4_a./repmat(trial_triplet_cat4_a(1,:), 3, 1);
errorbar(1:3, nanmean(trial_triplet_cat4_a, 2), sqrt(nanvar(trial_triplet_cat4_a, 0, 2)./sum(~isnan(trial_triplet_cat4_a),2)), 'r.-')
axis([0 4 0 35])

[h,p4a,stat4a,stat24a] = ttest(nanmean(trial_triplet_cat4_a([1,3], :),1) - trial_triplet_cat4_a(2,:));
title('Cat4_a')

subplot(2,2,2); hold on
% trial_triplet_cat6_a = trial_triplet_cat6_a./repmat(trial_triplet_cat6_a(1,:,:), 3,1,1);
%exclude subject 9: bc they are outlier
for i_sub_x = 1:length(additional_subjs_exclude)
    trial_triplet_cat6_a(:,:,additional_subjs_exclude(i_sub_x)) = nan;
end
for i_sub = 1:size(trial_triplet_cat6_a,3)
    compl_data = sum(sum(~isnan(trial_triplet_cat6_a(:,:,i_sub)),2)>0);
    if compl_data >= 3
        plot(1:3, nanmean(trial_triplet_cat6_a(:,:, i_sub), 2), '.-', 'Color', [.5 .5 .5]);
    end
end
trial_triplet_cat6_a = reshape(trial_triplet_cat6_a, 3, size(trial_triplet_cat6_a,2)*size(trial_triplet_cat6_a,3));
% trial_triplet_cat6_a = trial_triplet_cat6_a./repmat(trial_triplet_cat6_a(1,:), 3, 1);
errorbar((1:3), nanmean(trial_triplet_cat6_a, 2), sqrt(nanvar(trial_triplet_cat6_a, 0, 2)./sum(~isnan(trial_triplet_cat6_a),2)), 'g.-')
axis([0 4 0 35])

[h,p6a,stat6a,stat26a] = ttest(nanmean(trial_triplet_cat6_a([1,3], :),1) - trial_triplet_cat6_a(2,:));
title('Cat6_a')
%%

n_cat4 = nan(1, size(category_data,2));
for i_sub = 1:size(category_data,2)
    n_cat4(i_sub) = sum(category_data(trials_data(:, i_sub) > trials_b(1) & trials_data(:, i_sub) <= trials_b(2), i_sub) == 4);
end
n_cat6 = nan(1, size(category_data,2));
for i_sub = 1:size(category_data,2)
    cond1 = trials_data(:, i_sub) > trials_b(1) & trials_data(:, i_sub) <= trials_b(2);
    n_cat6(i_sub) = sum(category_data(cond1, i_sub) == 2 & probe_data(cond1, i_sub));
end
n_cat4_mode = mode(n_cat4);
n_cat6_mode = mode(n_cat6)-1;

trial_triplet_cat4_b = nan(3, n_cat4_mode, size(category_data, 2));
trial_triplet_cat6_b = nan(3, n_cat6_mode, size(category_data, 2));

for i_sub = 1:size(category_data,2)
    subset1 = trials_data(:, i_sub) > trials_b(1) & trials_data(:, i_sub) <= trials_b(2);
    reach_subset = reach_data(subset1, i_sub);
    reach_subset(reach_subset > 50 | reach_subset < -50) = nan;
    cat_subset = category_data(subset1, i_sub);
    probe_subset = probe_data(subset1, i_sub);
    cat2_inds = find(cat_subset == 2 & probe_subset == 0);
    cat4_inds = find(cat_subset == 4);
    cat6_inds = find(cat_subset == 2 & probe_subset == 1);
    cat2_p_inds = find(cat_subset == 2 & probe_off_subset == 1); % indicies of LPT trials that directly follow Aftereffect trials
    for i_4 = 1:size(trial_triplet_cat4_b, 2)
        trial_triplet_cat4_b(2, i_4, i_sub) = reach_subset(cat4_inds(i_4));
        k_before = max(cat2_inds(cat2_inds < cat4_inds(i_4)));
        k_after = min(cat2_inds(cat2_inds > cat4_inds(i_4)));
        if ~isempty(cat6_inds)
            if ismember(k_before, cat2_p_inds)
                % k_before is a post-aftereffect trial, omit it
                k_before = nan;
            end
            if ismember(k_after, cat2_p_inds)
                % k_after is a post-aftereffect trial, omit it
                k_after = nan;
            end
        end
        try
            if ~isnan(k_before)
                trial_triplet_cat4_b(1, i_4, i_sub) = reach_subset(k_before);
            end
            if ~isnan(k_after)
                trial_triplet_cat4_b(3, i_4, i_sub) = reach_subset(k_after);
            end
        catch
        end
    end
    if ~isempty(cat6_inds)
        for i_6 = 1:size(trial_triplet_cat6_b, 2)
            trial_triplet_cat6_b(2, i_6, i_sub) = reach_subset(cat6_inds(i_6));
            k_before = max(cat2_inds(cat2_inds < cat6_inds(i_6)));
            cat2_inds_after = cat2_inds(cat2_inds > cat6_inds(i_6));
            [k_after, inds_k_after] = min(cat2_inds_after);
            trial_triplet_cat6_b(1, i_6, i_sub) = reach_subset(k_before);
            try
                trial_triplet_cat6_b(3, i_6, i_sub) = reach_subset(cat2_inds_after(2));
            catch
            end
        end
    end
end

trial_triplet_cat4_b = trial_triplet_cat4_b(:, 2:end, :);
subplot(2,2,3); hold on
% trial_triplet_cat4_b = trial_triplet_cat4_b./repmat(trial_triplet_cat4_b(1,:,:), 3,1,1);
%exclude subject 9: bc they are outlier
for i_sub_x = 1:length(additional_subjs_exclude)
    trial_triplet_cat4_b(:,:,additional_subjs_exclude(i_sub_x)) = nan;
end
for i_sub = 1:size(trial_triplet_cat4_b,3)
    compl_data = sum(sum(~isnan(trial_triplet_cat4_b(:,:,i_sub)),2)>0);
    if compl_data >= 3
        plot(1:3, nanmean(trial_triplet_cat4_b(:,:, i_sub), 2), '.-', 'Color', [.5 .5 .5]);
    end
end
trial_triplet_cat4_b = reshape(trial_triplet_cat4_b, 3, size(trial_triplet_cat4_b,2)*size(trial_triplet_cat4_b,3));
% trial_triplet_cat4_a = trial_triplet_cat4_a./repmat(trial_triplet_cat4_a(1,:), 3, 1);
errorbar(1:3, nanmean(trial_triplet_cat4_b, 2), sqrt(nanvar(trial_triplet_cat4_b, 0, 2)./sum(~isnan(trial_triplet_cat4_b),2)), 'r.-')
axis([0 4 0 35])
[h,p4b,stat4b,stat24b] = ttest(nanmean(trial_triplet_cat4_b([1,3], :),1) - trial_triplet_cat4_b(2,:));
title('Cat4_b');

subplot(2,2,4); hold on
% trial_triplet_cat6_b = trial_triplet_cat6_b./repmat(trial_triplet_cat6_b(1,:,:), 3,1,1);
%exclude subject 9: bc they are outlier
for i_sub_x = 1:length(additional_subjs_exclude)
    trial_triplet_cat6_b(:,:,additional_subjs_exclude(i_sub_x)) = nan;
end
for i_sub = 1:size(trial_triplet_cat6_b,3)
    compl_data = sum(sum(~isnan(trial_triplet_cat6_b(:,:,i_sub)),2)>0);
    if compl_data >= 3
        plot(1:3, nanmean(trial_triplet_cat6_b(:,:, i_sub), 2), '.-', 'Color', [.5 .5 .5]);
    end
end
trial_triplet_cat6_b = reshape(trial_triplet_cat6_b, 3, size(trial_triplet_cat6_b,2)*size(trial_triplet_cat6_b,3));
% trial_triplet_cat6_b = trial_triplet_cat6_b./repmat(trial_triplet_cat6_b(1,:), 3, 1);
errorbar((1:3), nanmean(trial_triplet_cat6_b, 2), sqrt(nanvar(trial_triplet_cat6_b, 0, 2)./sum(~isnan(trial_triplet_cat6_b),2)), 'g.-')
axis([0 4 0 35])
[h,p6b,stat6b,stat26b] = ttest(nanmean(trial_triplet_cat6_b([1,3], :),1) - trial_triplet_cat6_b(2,:));
title('Cat6_b');
%% 
trial_triplet = {trial_triplet_cat4_a, trial_triplet_cat4_b, trial_triplet_cat6_a, trial_triplet_cat6_b};
stats = {p4a, p4b, p6a, p6b; stat24a, stat24b, stat26a, stat26b};