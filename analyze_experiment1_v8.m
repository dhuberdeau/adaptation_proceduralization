%% prep list of statistical figures:
stats_master_list = {};

%% Load data saved from manual trial exclusions

load('data_E1_ptFix_2TH_2cm_vManRev.mat')

subject_exclusions = [10 17]; % >50% of SPT trials discarded
all_subs = 1:size(category_data,2);
keep_subs = setdiff(all_subs, subject_exclusions);
category_temp = category_data;
pt_temp = pt_data;
reach_temp = reach_data;
trials_temp = trials_data;
probe_temp = probe_data;
clear category_data pt_data reach_data trials_data

category_data = category_temp(:, keep_subs);
pt_data = pt_temp(:, keep_subs);
reach_data = reach_temp(:, keep_subs);
trials_data = trials_temp(:, keep_subs);
probe_data = probe_temp(:, keep_subs);

N_SUBS = length(keep_subs);

%% define data by trial types
% cat* = "category *", where
%   cat1 = Long-PT trials to toward the target without the rotation
%   cat2 = Long-PT trials to toward the target with the rotation
%   cat3 = Short-PT trials to toward the target without the rotation
%   cat4 = Short-PT trials to toward the target with the rotation   

DISCARD_TH = 60; %throwaway outlier reaches

reach_data_cat2 = nan(1045, size(category_data,2));
reach_data_cat4 = nan(55, size(category_data,2));
reach_data_cat6 = nan(45, size(category_data,2));
for i_sub = 1:size(category_data,2)
    reach_data_cat2(:, i_sub) = reach_data(category_data(:, i_sub) == 2 & probe_data(:, i_sub) < .5, i_sub);
    reach_data_cat2(reach_data_cat2(:, i_sub) > DISCARD_TH | reach_data_cat2(:, i_sub) < -DISCARD_TH, i_sub) = nan;
    reach_data_cat4(:, i_sub) = reach_data(category_data(:, i_sub) == 4, i_sub);
    reach_data_cat4(reach_data_cat4(:, i_sub) > DISCARD_TH | reach_data_cat4(:, i_sub) < -DISCARD_TH, i_sub) = nan;
    reach_data_cat6(:, i_sub) = reach_data(category_data(:, i_sub) == 2 & probe_data(:, i_sub) > .5, i_sub);
    reach_data_cat6(reach_data_cat6(:, i_sub) > DISCARD_TH | reach_data_cat6(:, i_sub) < -DISCARD_TH, i_sub) = nan;
end

rotationA_onset_inds_cat2 = [100 240 390 540 690 840 985]; %rotation A onset
rotationB_onset_inds_cat2 = [140 290 440 590 740 890]; %rotation B onset
rotationN_onset_inds_cat2 = [190 340 490 640 790 940 1025]; %null onset
rotationA_onset_inds_cat4 = [25, 40];
rotationA_onset_inds_cat6 = [15, 30];

%% analyze rotation onset response

LOW_PT_TR_CNT = 2;

if LOW_PT_TR_CNT == 3
    rot_resp_wind_cat2 = [2 3 4 6 7 8 10 11 12]; %trials to include in measure of "rate" of adaptation for category 2 trials (excludes post-aftereffect (cat 6) trials)
elseif LOW_PT_TR_CNT == 2
    rot_resp_wind_cat2 = [2 3 4 6 7 8]; %trials to include in measure of "rate" of adaptation for category 2 trials (excludes post-aftereffect (cat 6) trials)
else
    error('No LOW_PT_TR_CNT valid');
end
rot_resp_wind_cat4_6 = 1:LOW_PT_TR_CNT; %trials to include in measure of "rate" of adaptation for category 4 trials

rotationA_onset_response_cat2 = nan(length(rotationA_onset_inds_cat2), size(reach_data_cat2,2));
rotationB_onset_response_cat2 = nan(length(rotationB_onset_inds_cat2), size(reach_data_cat2,2));
rotation_onset_response_cat4 = nan(length(rotationA_onset_inds_cat4), size(reach_data_cat4,2));
rotation_onset_response_cat6 = nan(length(rotationA_onset_inds_cat6), size(reach_data_cat6,2));

for i_rot = 1:length(rotationA_onset_inds_cat2)
    this_response = reach_data_cat2(rotationA_onset_inds_cat2(i_rot) + rot_resp_wind_cat2, :);
    for i_sub = 1:size(reach_data_cat2,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationA_onset_response_cat2(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationB_onset_inds_cat2)
    this_response = reach_data_cat2(rotationB_onset_inds_cat2(i_rot) + rot_resp_wind_cat2, :);
    for i_sub = 1:size(reach_data_cat2,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationB_onset_response_cat2(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationA_onset_inds_cat4)
    this_response = reach_data_cat4(rotationA_onset_inds_cat4(i_rot) + rot_resp_wind_cat4_6, :);
    for i_sub = 1:size(reach_data_cat4,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotation_onset_response_cat4(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationA_onset_inds_cat6)
    this_response = reach_data_cat6(rotationA_onset_inds_cat6(i_rot) + rot_resp_wind_cat4_6, :);
    for i_sub = 1:size(reach_data_cat6,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotation_onset_response_cat6(i_rot,:) = nanmean(this_response, 1);
end

%% Do anova and multiple comparisons (rotation onset response)
onset_response_matrix = [rotationA_onset_response_cat2(1,:)', rotation_onset_response_cat4(1,:)', rotation_onset_response_cat6(1,:)'];
for i_sub = 1:size(onset_response_matrix, 1)
    if sum(isnan(onset_response_matrix(i_sub, :))) > 0
        onset_response_matrix(i_sub, :) = nan(1,3);
    end
end
[a,b,c] = anova1(onset_response_matrix);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate anova p';
stats_master_list{size(stats_master_list,1), 2} = a;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate anova F';
stats_master_list{size(stats_master_list,1), 2} = b{2,5};

[h_2_4, p_2_4, ~, s_2_4] = ttest(onset_response_matrix(:,1) - onset_response_matrix(:,2));
[h_2_6, p_2_6, ~, s_2_6] = ttest(onset_response_matrix(:,1) - onset_response_matrix(:,3));
[h_4_6, p_4_6, ~, s_4_6] = ttest(onset_response_matrix(:,2) - onset_response_matrix(:,3));

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate cat2_vs_cat4 p';
stats_master_list{size(stats_master_list,1), 2} = p_2_4;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate cat2_vs_cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = s_2_4.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate cat2_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = p_2_6;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate cat2_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = s_2_6.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate cat4_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = p_4_6;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 rate cat4_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = s_4_6.tstat;
%% analyze asymptote:
if LOW_PT_TR_CNT == 3
    asym_wind_cat2 = [30 31 32 34 35 36 38 39 40]; %trials to include in measure of "rate" of adaptation for category 2 trials (excludes post-aftereffect (cat 6) trials)
elseif LOW_PT_TR_CNT == 2
    asym_wind_cat2 = [34 35 36 38 39 40]; %trials to include in measure of "rate" of adaptation for category 2 trials (excludes post-aftereffect (cat 6) trials)
else
    error('No LOW_PT_TR_CNT valid');
end
% asym_wind_cat2 = [28 30 31 32 34 35 36 38 39 40];
asym_wind_cat4_6 = (1:LOW_PT_TR_CNT)+8 - 1;

rotationA_asym_cat2 = nan(length(rotationA_onset_inds_cat2), size(reach_data_cat2,2));
rotationB_asym_cat2 = nan(length(rotationB_onset_inds_cat2), size(reach_data_cat2,2));
rotation_asym_cat4 = nan(length(rotationA_onset_inds_cat4), size(reach_data_cat4,2));
rotation_asym_cat6 = nan(length(rotationA_onset_inds_cat6), size(reach_data_cat6,2));

for i_rot = 1:length(rotationA_onset_inds_cat2)
    this_response = reach_data_cat2(rotationA_onset_inds_cat2(i_rot) + asym_wind_cat2, :);
    for i_sub = 1:size(reach_data_cat2,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationA_asym_cat2(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationB_onset_inds_cat2)
    this_response = reach_data_cat2(rotationB_onset_inds_cat2(i_rot) + asym_wind_cat2, :);
    for i_sub = 1:size(reach_data_cat2,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationB_asym_cat2(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationA_onset_inds_cat4)
    this_response = reach_data_cat4(rotationA_onset_inds_cat4(i_rot) + asym_wind_cat4_6, :);
    for i_sub = 1:size(reach_data_cat4,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotation_asym_cat4(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationA_onset_inds_cat6)
    this_response = reach_data_cat6(rotationA_onset_inds_cat6(i_rot) + asym_wind_cat4_6, :);
    for i_sub = 1:size(reach_data_cat6,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotation_asym_cat6(i_rot,:) = nanmean(this_response, 1);
end

%% Do anova and multiple comparisons (asymptote)
asym_matrix = [rotationA_asym_cat2(1,:)', rotation_asym_cat4(1,:)', rotation_asym_cat6(1,:)'];
for i_sub = 1:size(asym_matrix, 1)
    if sum(isnan(asym_matrix(i_sub, :))) > 0
        asym_matrix(i_sub, :) = nan(1,3);
    end
end
[a_a,b_a,c_a] = anova1(asym_matrix);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym anova p';
stats_master_list{size(stats_master_list,1), 2} = a_a;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym anova F';
stats_master_list{size(stats_master_list,1), 2} = b_a{2,5};

[h_2_4_a, p_2_4_a, ~, s_2_4_a] = ttest(asym_matrix(:,1) - asym_matrix(:,2));
[h_2_6_a, p_2_6_a, ~, s_2_6_a] = ttest(asym_matrix(:,1) - asym_matrix(:,3));
[h_4_6_a, p_4_6_a, ~, s_4_6_a] = ttest(asym_matrix(:,2) - asym_matrix(:,3));

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym cat2_vs_cat4 p';
stats_master_list{size(stats_master_list,1), 2} = p_2_4_a;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym cat2_vs_cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = s_2_4_a.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym cat2_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = p_2_6_a;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym cat2_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = s_2_6_a.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym cat4_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = p_4_6_a;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 asym cat4_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = s_4_6_a.tstat;

%% nearest neighbor analysis:

trials1_a = [280, 316]; % trial numbers - rate, cycle 1
trials1_b = [364 400];% trial numbers - asym, cycle 1
[~, stats1] = nearest_neighbor_analysis_fnc(category_data, reach_data, trials_data, probe_data, trials1_a, trials1_b, subject_exclusions);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborRate cat2_vs_cat4 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,1};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborRate cat2_vs_cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,1}.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborAsym cat2_vs_cat4 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,2};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborAsym cat2_vs_cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,2}.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborRate cat2_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,3};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborRate cat2_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,3}.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborAsym cat2_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,4};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 nearestNeighborAsym cat2_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,4}.tstat;


trials1_a = [2110, 2146]; % trial numbers - rate, cycle 7
trials1_b = [2194 2230]; % trial numbers - asym, cycle 7
[~, stats1] = nearest_neighbor_analysis_fnc(category_data, reach_data, trials_data, probe_data, trials1_a, trials1_b, subject_exclusions);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborRate cat2_vs_cat4 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,1};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborRate cat2_vs_cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,1}.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborAsym cat2_vs_cat4 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,2};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborAsym cat2_vs_cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,2}.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborRate cat2_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,3};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborRate cat2_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,3}.tstat;

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborAsym cat2_vs_cat6 p';
stats_master_list{size(stats_master_list,1), 2} = stats1{1,4};
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 nearestNeighborAsym cat2_vs_cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = stats1{2,4}.tstat;

%% post-aftereffect regression trials
post_cat6_trials = [5 9 13 17 21 25 29 33 37];
pre_cat6_trials = post_cat6_trials - 1;

post_tr = nan(length(post_cat6_trials), size(reach_data_cat2,2));
pre_tr = nan(length(pre_cat6_trials), size(reach_data_cat2,2));
for i_tr = 1:length(post_cat6_trials)
    post_tr(i_tr, :) = reach_data_cat2(rotationA_onset_inds_cat2(1) + post_cat6_trials(i_tr), :);
    pre_tr(i_tr, :) = reach_data_cat2(rotationA_onset_inds_cat2(1) + pre_cat6_trials(i_tr), :);
end
mean_diff_regression_cyc1 = nanmean(post_tr - pre_tr,1);
[~,p_AF_regression_1,~,stats_regression_1] = ttest(mean_diff_regression_cyc1);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 all regression p';
stats_master_list{size(stats_master_list,1), 2} = p_AF_regression_1;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot1 all regression tstat';
stats_master_list{size(stats_master_list,1), 2} = stats_regression_1.tstat;

post_cat6_trials = [5 9 13 17 21 25 29 33 37];
pre_cat6_trials = post_cat6_trials - 1;

post_tr = nan(length(post_cat6_trials), size(reach_data_cat2,2));
pre_tr = nan(length(pre_cat6_trials), size(reach_data_cat2,2));
for i_tr = 1:length(post_cat6_trials)
    post_tr(i_tr, :) = reach_data_cat2(rotationA_onset_inds_cat2(7) + post_cat6_trials(i_tr), :);
    pre_tr(i_tr, :) = reach_data_cat2(rotationA_onset_inds_cat2(7) + pre_cat6_trials(i_tr), :);
end
mean_diff_regression_cyc7 = nanmean(post_tr - pre_tr,1);
[~,p_AF_regression_7,~,stats_regression_7] = ttest(mean_diff_regression_cyc7);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 all regression p';
stats_master_list{size(stats_master_list,1), 2} = p_AF_regression_7;
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rot7 all regression tstat';
stats_master_list{size(stats_master_list,1), 2} = stats_regression_7.tstat;

%% analyze savings interaction in rate

onset_response_matrix_a = [rotationA_onset_response_cat2(1,:)', rotation_onset_response_cat4(1,:)', rotation_onset_response_cat6(1,:)'];
onset_response_matrix_b = [rotationA_onset_response_cat2(end,:)', rotation_onset_response_cat4(end,:)', rotation_onset_response_cat6(end,:)'];
nan_subs_a = nan(1,N_SUBS); n_nan_subs = 0;
for i_sub = 1:(size(onset_response_matrix_a, 1))
    if sum(isnan(onset_response_matrix_a(i_sub, :))) > 0
        onset_response_matrix_a(i_sub, :) = nan(1,3);
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_a(n_nan_subs) = i_sub;
    end
end
nan_subs_b = nan(1,N_SUBS); n_nan_subs = 0;
for i_sub = 1:(size(onset_response_matrix_b, 1))
    if sum(isnan(onset_response_matrix_b(i_sub, :))) > 0
        onset_response_matrix_b(i_sub, :) = nan(1,3);
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_b(n_nan_subs) = i_sub;
    end
end

valid_subs = setdiff(1:N_SUBS, union(nan_subs_a, nan_subs_b));
onset_response_matrix_anova_valid = [onset_response_matrix_a(valid_subs,:); onset_response_matrix_b(valid_subs,:)];

[a_sav, b_sav, c_sav] = anova2(onset_response_matrix_anova_valid, 2);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate all interaction p';
stats_master_list{size(stats_master_list,1), 2} = a_sav(3);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate all interaction F';
stats_master_list{size(stats_master_list,1), 2} = b_sav{4, 5};
%% analyze savings interaction in asymptote
asym_matrix_a = [rotationA_asym_cat2(1,:)', rotation_asym_cat4(1,:)', rotation_asym_cat6(1,:)'];
asym_matrix_b = [rotationA_asym_cat2(end,:)', rotation_asym_cat4(end,:)', rotation_asym_cat6(end,:)'];
nan_subs_a = nan(1,N_SUBS); n_nan_subs = 0;
for i_sub = 1:(size(asym_matrix_a, 1))
    if sum(isnan(asym_matrix_a(i_sub, :))) > 0
        asym_matrix_a(i_sub, :) = nan(1,3);
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_a(n_nan_subs) = i_sub;
    end
end
nan_subs_b = nan(1,N_SUBS); n_nan_subs = 0;
for i_sub = 1:(size(asym_matrix_b, 1))
    if sum(isnan(asym_matrix_b(i_sub, :))) > 0
        asym_matrix_b(i_sub, :) = nan(1,3);
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_b(n_nan_subs) = i_sub;
    end
end

valid_subs = setdiff(1:N_SUBS, union(nan_subs_a, nan_subs_b));
asym_matrix_anova_valid = [asym_matrix_a(valid_subs,:); asym_matrix_b(valid_subs,:)];
% asym_matrix_anova_valid = [asym_matrix_a; asym_matrix_b];

[a_sav, b_sav, c_sav] = anova2(asym_matrix_anova_valid, 2);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym all interaction p';
stats_master_list{size(stats_master_list,1), 2} = a_sav(3);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym all interaction F';
stats_master_list{size(stats_master_list,1), 2} = b_sav{4, 5};

%% analyze savings in rate and asymptote
savings_matrix_rotation_response = [rotationA_onset_response_cat2(1,:)' - rotationA_onset_response_cat2(end,:)',...
    rotation_onset_response_cat4(1,:)' - rotation_onset_response_cat4(end,:)',...
    rotation_onset_response_cat6(1,:)' - rotation_onset_response_cat6(end,:)'];

[h_savings, p_savings, ~, d_savings] = ttest(savings_matrix_rotation_response);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate cat2 p';
stats_master_list{size(stats_master_list,1), 2} = p_savings(1);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate cat2 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_savings.tstat(1);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate cat4 p';
stats_master_list{size(stats_master_list,1), 2} = p_savings(2);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_savings.tstat(2);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate cat6 p';
stats_master_list{size(stats_master_list,1), 2} = p_savings(3);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings rate cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_savings.tstat(3);


savings_matrix_asym = [rotationA_asym_cat2(1,:)' - rotationA_asym_cat2(end,:)',...
    rotation_asym_cat4(1,:)' - rotation_asym_cat4(end,:)',...
    rotation_asym_cat6(1,:)' - rotation_asym_cat6(end,:)'];

[h_asym, p_asym, ~, d_asym] = ttest(savings_matrix_asym);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym cat2 p';
stats_master_list{size(stats_master_list,1), 2} = p_asym(1);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym cat2 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_asym.tstat(1);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym cat4 p';
stats_master_list{size(stats_master_list,1), 2} = p_asym(2);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym cat4 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_asym.tstat(2);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym cat6 p';
stats_master_list{size(stats_master_list,1), 2} = p_asym(3);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'savings asym cat6 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_asym.tstat(3);


savings_matrix_rotation_response = rotationB_onset_response_cat2(1,:)' - rotationB_onset_response_cat2(end,:)';
[h_savings_B, p_savings_B, ~, d_savings_B] = ttest(savings_matrix_rotation_response);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rotation B savings rate cat2 p';
stats_master_list{size(stats_master_list,1), 2} = p_savings_B(1);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rotation B savings rate cat2 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_savings_B.tstat(1);

savings_matrix_rotation_response = rotationB_asym_cat2(1,:)' - rotationB_asym_cat2(end,:)';
[h_savings_aB, p_savings_aB, ~, d_savings_aB] = ttest(savings_matrix_rotation_response);

stats_master_list{size(stats_master_list,1) + 1, 1} = 'rotation B savings asym cat2 p';
stats_master_list{size(stats_master_list,1), 2} = p_savings_aB(1);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'rotation B savings asym cat2 tstat';
stats_master_list{size(stats_master_list,1), 2} = d_savings_aB.tstat(1);



%% pre-rotation onset

rot_resp_wind_cat2 = 1:8;
rot_resp_wind_cat4_6 = 1:LOW_PT_TR_CNT;
asym_wind_cat2 = 33:40;
asym_wind_cat4_6 = (1:LOW_PT_TR_CNT)+8 - 1;

rotationA_pre_response_cat2 = nan(length(rotationA_onset_inds_cat2), size(reach_data_cat2,2));
rotationB_pre_response_cat2 = nan(length(rotationB_onset_inds_cat2), size(reach_data_cat2,2));
rotation_pre_response_cat4 = nan(length(rotationA_onset_inds_cat4), size(reach_data_cat4,2));
rotation_pre_response_cat6 = nan(length(rotationA_onset_inds_cat6), size(reach_data_cat6,2));

for i_rot = 1:length(rotationA_onset_inds_cat2)
    rotationA_pre_response_cat2(i_rot,:) = nanmean(reach_data_cat2(rotationA_onset_inds_cat2(i_rot) - rot_resp_wind_cat2, :), 1);
end

for i_rot = 1:length(rotationB_onset_inds_cat2)
    rotationB_pre_response_cat2(i_rot,:) = nanmean(reach_data_cat2(rotationB_onset_inds_cat2(i_rot) - rot_resp_wind_cat2, :), 1);
end

for i_rot = 1:length(rotationA_onset_inds_cat4)
    rotation_pre_response_cat4(i_rot,:) = nanmean(reach_data_cat4(rotationA_onset_inds_cat4(i_rot) - rot_resp_wind_cat4_6, :), 1);
end

for i_rot = 1:length(rotationA_onset_inds_cat6)
    rotation_pre_response_cat6(i_rot,:) = nanmean(reach_data_cat6(rotationA_onset_inds_cat6(i_rot) - rot_resp_wind_cat4_6, :), 1);
end

%% multiple comparisons of pre-rotation period
[a2,b2,c2] = anova1([rotationA_pre_response_cat2(1,:)', rotation_pre_response_cat4(1,:)', rotation_pre_response_cat6(1,:)'])
[h_re,p_re,j_junk,s_re] = ttest(rotationA_pre_response_cat2(7,:) - rotationA_pre_response_cat2(1,:))
stats_master_list{size(stats_master_list,1) + 1, 1} = 'pre-rotation cycles1_vs_7 cat2 p';
stats_master_list{size(stats_master_list,1), 2} = p_re(1);
stats_master_list{size(stats_master_list,1) + 1, 1} = 'pre-rotation cycles1_vs_7 cat2 t-test';
stats_master_list{size(stats_master_list,1), 2} = s_re.tstat(1);