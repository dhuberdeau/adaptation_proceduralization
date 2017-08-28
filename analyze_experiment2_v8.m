%% prep list of statistical figures:
stats_master_list_2 = {};

%% Load data saved from manual trial exclusions

load('data_E2_ptFix_2TH_2cm_vManRev_fixed.mat')

subject_exclusions = [7 9 17]; %bc >half of short-PT trials are faulty
all_subs = 1:size(category_data,2);
keep_subs = setdiff(all_subs, subject_exclusions);
category_temp = category_data;
pt_temp = pt_data;
reach_temp = reach_data;
trials_temp = trials_data;
clear category_data pt_data reach_data trials_data

category_data = category_temp(:, keep_subs);
pt_data = pt_temp(:, keep_subs);
reach_data = reach_temp(:, keep_subs);
trials_data = trials_temp(:, keep_subs);

%% define data by trial types
% cat* = "category *", where
%   cat1 = Long-PT trials to toward the target without the rotation
%   cat2 = Long-PT trials to toward the target with the rotation
%   cat3 = Short-PT trials to toward the target without the rotation
%   cat4 = Short-PT trials to toward the target with the rotation   

DISCARD_TH = 60; %throwaway outlier reaches

N_CAT2_TRS = 880;
N_CAT4_TRS = 220;

reach_data_cat2 = nan(N_CAT2_TRS, size(category_data,2));
reach_data_cat4 = nan(N_CAT4_TRS, size(category_data,2));
for i_sub = 1:size(category_data,2)
    reach_data_cat2(:, i_sub) = reach_data(category_data(:, i_sub) == 2, i_sub);
    reach_data_cat2(reach_data_cat2(:, i_sub) > DISCARD_TH | reach_data_cat2(:, i_sub) < -DISCARD_TH, i_sub) = nan;
    reach_data_cat4(:, i_sub) = reach_data(category_data(:, i_sub) == 4, i_sub);
    reach_data_cat4(reach_data_cat4(:, i_sub) > DISCARD_TH | reach_data_cat4(:, i_sub) < -DISCARD_TH, i_sub) = nan;
end

rotation_change_inds = 100:40:N_CAT2_TRS; 

rotationA_onset_inds_cat2 = rotation_change_inds(1:3:20); %rotation A onset
rotationB_onset_inds_cat2 = rotation_change_inds(2:3:17); %rotation B onset
rotationN_onset_inds_cat2 = rotation_change_inds(3:3:20); %rotation null onset

rotation_change_inds_cat4 = 25:10:N_CAT4_TRS;
rotationA_onset_inds_cat4 = rotation_change_inds_cat4(1:3:20); %rotation A onset
rotationB_onset_inds_cat4 = rotation_change_inds_cat4(2:3:17); %rotation B onset
rotationN_onset_inds_cat4 = rotation_change_inds_cat4(3:3:20); %rotation null onset

%% analyze rotation onset response

LOW_PT_TR_CNT = 2;
rot_resp_wind_cat2 = 2:8; %trials to include in measure of "rate" of adaptation for category 2 trials (excludes post-aftereffect (cat 6) trials)

rot_resp_wind_cat4 = 1:LOW_PT_TR_CNT; %trials to include in measure of "rate" of adaptation for category 4 trials

rotationA_onset_response_cat2 = nan(length(rotationA_onset_inds_cat2), size(reach_data_cat2,2));
rotationB_onset_response_cat2 = nan(length(rotationB_onset_inds_cat2), size(reach_data_cat2,2));
rotationA_onset_response_cat4 = nan(length(rotationA_onset_inds_cat4), size(reach_data_cat4,2));
rotationB_onset_response_cat4 = nan(length(rotationB_onset_inds_cat4), size(reach_data_cat4,2));


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
    this_response = reach_data_cat4(rotationA_onset_inds_cat4(i_rot) + rot_resp_wind_cat4, :);
    for i_sub = 1:size(reach_data_cat4,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationA_onset_response_cat4(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationB_onset_inds_cat4)
    this_response = reach_data_cat4(rotationB_onset_inds_cat4(i_rot) + rot_resp_wind_cat4, :);
    for i_sub = 1:size(reach_data_cat4,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationB_onset_response_cat4(i_rot,:) = nanmean(this_response, 1);
end

%% Do anova and multiple comparisons (rotation onset response)
onset_response_matrix = [rotationA_onset_response_cat2(1,:)', rotationA_onset_response_cat4(1,:)'];
for i_sub = 1:size(onset_response_matrix, 1)
    if sum(isnan(onset_response_matrix(i_sub, :))) > 0
        onset_response_matrix(i_sub, :) = nan(1,size(onset_response_matrix,2));
    end
end

[h_2_4, p_2_4, ~, s_2_4] = ttest(onset_response_matrix(:,1) - onset_response_matrix(:,2));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 A rate cat2_vs_cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_2_4;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 A rate cat2_vs_cat4 tstat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_2_4.tstat;

% rotation B
onset_response_matrix = [rotationB_onset_response_cat2(1,:)', rotationB_onset_response_cat4(1,:)'];
for i_sub = 1:size(onset_response_matrix, 1)
    if sum(isnan(onset_response_matrix(i_sub, :))) > 0
        onset_response_matrix(i_sub, :) = nan(1,size(onset_response_matrix,2));
    end
end

[h_2_4, p_2_4, ~, s_2_4] = ttest(onset_response_matrix(:,1) - onset_response_matrix(:,2));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 B rate cat2_vs_cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_2_4;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 B rate cat2_vs_cat4 tstat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_2_4.tstat;
%% analyze asymptote:

asym_wind_cat2 = 34:40; %trials to include in measure of "rate" of adaptation for category 2 trials
asym_wind_cat4 = (1:LOW_PT_TR_CNT) + 8 - 1;

rotationA_asym_cat2 = nan(length(rotationA_onset_inds_cat2), size(reach_data_cat2,2));
rotationB_asym_cat2 = nan(length(rotationB_onset_inds_cat2), size(reach_data_cat2,2));
rotationA_asym_cat4 = nan(length(rotationA_onset_inds_cat4), size(reach_data_cat4,2));
rotationB_asym_cat4 = nan(length(rotationB_onset_inds_cat4), size(reach_data_cat4,2));

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
    this_response = reach_data_cat4(rotationA_onset_inds_cat4(i_rot) + asym_wind_cat4, :);
    for i_sub = 1:size(reach_data_cat4,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationA_asym_cat4(i_rot,:) = nanmean(this_response, 1);
end

for i_rot = 1:length(rotationB_onset_inds_cat4)
    this_response = reach_data_cat4(rotationB_onset_inds_cat4(i_rot) + asym_wind_cat4, :);
    for i_sub = 1:size(reach_data_cat4,2)
        if sum(isnan(this_response(:, i_sub))) > size(this_response,1)/2
            this_response(:, i_sub) = nan;
        end
    end
    rotationB_asym_cat4(i_rot,:) = nanmean(this_response, 1);
end

%% Do anova and multiple comparisons (asymptote)
asym_matrix = [rotationA_asym_cat2(1,:)', rotationA_asym_cat4(1,:)'];
for i_sub = 1:size(asym_matrix, 1)
    if sum(isnan(asym_matrix(i_sub, :))) > 0
        asym_matrix(i_sub, :) = nan(1,size(asym_matrix,2));
    end
end

[h_2_4_a, p_2_4_a, ~, s_2_4_a] = ttest(asym_matrix(:,1) - asym_matrix(:,2));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 A asym cat2_vs_cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_2_4_a;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 A asym cat2_vs_cat4 tstat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_2_4_a.tstat;

asym_matrix = [rotationB_asym_cat2(1,:)', rotationB_asym_cat4(1,:)'];
for i_sub = 1:size(asym_matrix, 1)
    if sum(isnan(asym_matrix(i_sub, :))) > 0
        asym_matrix(i_sub, :) = nan(1,size(asym_matrix,2));
    end
end

[h_2_4_a, p_2_4_a, ~, s_2_4_a] = ttest(asym_matrix(:,1) - asym_matrix(:,2));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 B asym cat2_vs_cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_2_4_a;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rot1 B asym cat2_vs_cat4 tstat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_2_4_a.tstat;


%% analyze savings interaction in rate, rotation A

rotationA_onset_response_cat2_purge = rotationA_onset_response_cat2([1, 7], :);
rotationA_onset_response_cat4_purge = rotationA_onset_response_cat4([1, 7], :);

nan_subs_cat2 = nan(1,size(rotationA_onset_response_cat2_purge,2)); n_nan_subs = 0;
for i_sub = 1:(size(rotationA_onset_response_cat2_purge, 2))
    if sum(isnan(rotationA_onset_response_cat2_purge(:,i_sub))) > 0
        rotationA_onset_response_cat2_purge(:,i_sub) = nan;
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_cat2(n_nan_subs) = i_sub;
    end
end
nan_subs_cat4 = nan(1,size(rotationA_onset_response_cat4_purge,2)); n_nan_subs = 0;
for i_sub = 1:(size(rotationA_onset_response_cat4_purge, 2))
    if sum(isnan(rotationA_onset_response_cat4_purge(:,i_sub))) > 0
        rotationA_onset_response_cat4_purge(:,i_sub) = nan;
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_cat4(n_nan_subs) = i_sub;
    end
end

valid_subs = setdiff(1:size(rotationA_onset_response_cat2,2), union(nan_subs_cat2, nan_subs_cat4));
onset_response_matrix_A_anova_valid = [rotationA_onset_response_cat2_purge(:, valid_subs)'; rotationA_onset_response_cat4_purge(:, valid_subs)'];

[a_sav, b_sav, c_sav] = anova2(onset_response_matrix_A_anova_valid, length(valid_subs));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotation savings rate all interaction p';
stats_master_list_2{size(stats_master_list_2,1), 2} = a_sav(3);
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'savings rate all anova F';
stats_master_list_2{size(stats_master_list_2,1), 2} = b_sav{4, 5};

% save data for analysis in R
mat1 = [rotationA_onset_response_cat2'; rotationA_onset_response_cat4']; %rotation rate data
mat2 = repmat(1:size(rotationA_onset_response_cat2,1), size(rotationA_onset_response_cat2,2) + size(rotationA_onset_response_cat4,2), 1); % cycle number factor
mat3 = [zeros(size(rotationA_onset_response_cat2')); ones(size(rotationA_onset_response_cat4'))]; % trial type factor
mat4 = repmat([(1:size(rotationA_onset_response_cat2,2))'; (1:size(rotationA_onset_response_cat4,2))'], 1, size(rotationA_onset_response_cat2,1)); % subject number
output_matrix_A = mat1(:); % response (rate data)
factor_matrix_A = [mat2(:), mat3(:), mat4(:)]; % factors: cycle, type, subject

csvwrite(['/Users/david/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_proceduralization/', 'rotationA_output_matrix'], output_matrix_A);
csvwrite(['/Users/david/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_proceduralization/', 'rotationA_factor_matrix'], factor_matrix_A);

%% analyze savings interaction in rate, rotation B
rotationB_onset_response_cat2_purge = rotationB_onset_response_cat2([1,6], :);
rotationB_onset_response_cat4_purge = rotationB_onset_response_cat4([1,6], :);

nan_subs_cat2 = nan(1,size(rotationB_onset_response_cat2_purge,2)); n_nan_subs = 0;
for i_sub = 1:(size(rotationB_onset_response_cat2_purge, 2))
    if sum(isnan(rotationB_onset_response_cat2_purge(:,i_sub))) > 0
        rotationB_onset_response_cat2_purge(:,i_sub) = nan;
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_cat2(n_nan_subs) = i_sub;
    end
end
nan_subs_cat4 = nan(1,size(rotationB_onset_response_cat4_purge,2)); n_nan_subs = 0;
for i_sub = 1:(size(rotationB_onset_response_cat4_purge, 2))
    if sum(isnan(rotationB_onset_response_cat4_purge(:,i_sub))) > 0
        rotationB_onset_response_cat4_purge(:,i_sub) = nan;
        n_nan_subs = 1 + n_nan_subs;
        nan_subs_cat4(n_nan_subs) = i_sub;
    end
end

valid_subs = setdiff(1:size(rotationB_onset_response_cat2,2), union(nan_subs_cat2, nan_subs_cat4));
onset_response_matrix_B_anova_valid = [rotationB_onset_response_cat2_purge(:, valid_subs)'; rotationB_onset_response_cat4_purge(:, valid_subs)'];

[a_sav, b_sav, c_sav] = anova2(onset_response_matrix_B_anova_valid, length(valid_subs));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'savings rate all interaction p';
stats_master_list_2{size(stats_master_list_2,1), 2} = a_sav(3);
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'savings rate all interaction F';
stats_master_list_2{size(stats_master_list_2,1), 2} = b_sav{4, 5};

% save data for analysis in R
mat1 = [rotationB_onset_response_cat2'; rotationB_onset_response_cat4']; %rotation rate data
mat2 = repmat(1:size(rotationB_onset_response_cat2,1), size(rotationB_onset_response_cat2,2) + size(rotationB_onset_response_cat4,2), 1); % cycle number factor
mat3 = [zeros(size(rotationB_onset_response_cat2')); ones(size(rotationB_onset_response_cat4'))]; % trial type factor
mat4 = repmat([(1:size(rotationB_onset_response_cat2,2))'; (1:size(rotationB_onset_response_cat4,2))'], 1, size(rotationB_onset_response_cat2,1)); % subject number
output_matrix_B = mat1(:); % response (rate data)
factor_matrix_B = [mat2(:), mat3(:), mat4(:)]; % factors: cycle, type, subject

csvwrite(['/Users/david/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_proceduralization/', 'rotationB_output_matrix'], output_matrix_B);
csvwrite(['/Users/david/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_proceduralization/', 'rotationB_factor_matrix'], factor_matrix_B);

%% cat2 vs. cat4 Cycle 6 rotation B
[h_cyc6, p_cyc6, ~, s_cyc6] = ttest(rotationB_onset_response_cat2(6,:) - rotationB_onset_response_cat4(6,:));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationB rate cycle6 cat2_vs_cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_cyc6;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationB rate cycle6 cat2_vs_cat4 t-stat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_cyc6.tstat;


%%  savings cyc1 vs 7 rotationA 
[~, p_sav, ~, s_sav] = ttest(rotationA_onset_response_cat2(1,:) - rotationA_onset_response_cat2(7,:));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationA savings cycle1_vs_7 cat2 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_sav;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationA savings cycle1_vs_7 cat2 t-stat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_sav.tstat;

[~, p_sav, ~, s_sav] = ttest(rotationA_onset_response_cat4(1,:) - rotationA_onset_response_cat4(7,:));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationA savings cycle1_vs_7 cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_sav;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationA savings cycle1_vs_7 cat4 t-stat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_sav.tstat;

%%  savings cyc1 vs 6 rotationB 
[~, p_sav, ~, s_sav] = ttest(rotationB_onset_response_cat2(1,:) - rotationB_onset_response_cat2(6,:));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationB rate cycle1_vs_6 cat2 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_sav;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationB rate cycle1_vs_6 cat2_vs_cat4 t-stat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_sav.tstat;

[~, p_sav, ~, s_sav] = ttest(rotationB_onset_response_cat4(1,:) - rotationB_onset_response_cat4(6,:));

stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationB rate cycle1_vs_6 cat4 p';
stats_master_list_2{size(stats_master_list_2,1), 2} = p_sav;
stats_master_list_2{size(stats_master_list_2,1) + 1, 1} = 'rotationB rate cycle1_vs_6 cat4 t-stat';
stats_master_list_2{size(stats_master_list_2,1), 2} = s_sav.tstat;

