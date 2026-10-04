clear all
close all
clc

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();


% /// Pipeline script #7: analyze hemispheric differences and compare groups /// 

%%  Set user-defined parameters

% Define your two groups
ctrl_type = 'naive'; 
exp_type  = 'rws';

% Base directory (common part)
base_root = paths.data;  

% Construct full paths
ctrl_dir = fullfile(base_root, ctrl_type);
exp_dir  = fullfile(base_root, exp_type);

%% Allen atlas setup

allenDir = paths.atlas;
addpath(allenDir);
AllenFile = fullfile(allenDir, 'annotation_10.nii.gz');
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2),:,:);
brainMask = AllenCrop > 0;
half_atlas = AllenCrop(:,:,1:end);

%% Load scaled_diff_4d_new for both groups

% Load control group data 
scaled_diff_4d_new_ctrl = load(fullfile(ctrl_dir, 'scaled_diff_4d_new.mat'), 'scaled_diff_4d_new');
scaled_diff_4d_new_ctrl = scaled_diff_4d_new_ctrl.scaled_diff_4d_new;

% Load experimental group data 
scaled_diff_4d_new_exp = load(fullfile(exp_dir, 'scaled_diff_4d_new.mat'), 'scaled_diff_4d_new');
scaled_diff_4d_new_exp = scaled_diff_4d_new_exp.scaled_diff_4d_new;

%% Compute LR differences and group differences

% Compute LR difference
[lr_diff_ctrl, lr_sum_ctrl] = compute_lr_stats(scaled_diff_4d_new_ctrl);
[lr_diff_exp,  lr_sum_exp]  = compute_lr_stats(scaled_diff_4d_new_exp);
avg_lr_diff_ctrl = nanmean(abs(lr_diff_ctrl), 4); %#ok<NANMEAN>
avg_lr_sum_ctrl  = nanmean(abs(lr_sum_ctrl),  4);
avg_lr_diff_exp  = nanmean(abs(lr_diff_exp),  4);
avg_lr_sum_exp   = nanmean(abs(lr_sum_exp),   4);
avg_lr_diff_groupdiff = avg_lr_diff_exp - avg_lr_diff_ctrl;
avg_lr_sum_groupdiff  = avg_lr_sum_exp  - avg_lr_sum_ctrl;

% Crop mask and atlas to actual data size
brainMask_cropped = brainMask(:, :, 1:size(avg_lr_diff_ctrl, 3));

% Generate average ctrl group video
write_lr_video(avg_lr_diff_ctrl, avg_lr_sum_ctrl, half_atlas, brainMask_cropped, ...
    ctrl_dir, ['lr_diff_sum_' ctrl_type '.mp4'], [0, 0.033], ctrl_type, ...
    'LR abs difference (diff) average', 'LR abs difference (sum) average');

% Generate average exp group video
write_lr_video(avg_lr_diff_exp, avg_lr_sum_exp, half_atlas, brainMask_cropped, ...
    exp_dir, ['lr_diff_sum_' exp_type '.mp4'], [0, 0.033], exp_type, ...
    'LR abs difference (diff) average', 'LR abs difference (sum) average');

% Generate average group difference video
write_lr_video(avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, brainMask_cropped, ...
    base_root, ['lr_diff_sum_groupdiff_' ctrl_type '_' exp_type '.mp4'], [-0.033, 0.033], ...
    [exp_type ' - ' ctrl_type], 'LR abs diff groupdiff', 'LR abs sum groupdiff');

%% Plot t-scored videos

sem_lr_diff_ctrl = nanstd(abs(lr_diff_ctrl), [], 4) ./ sqrt(sum(~isnan(lr_diff_ctrl),4));
sem_lr_sum_ctrl  = nanstd(abs(lr_sum_ctrl),  [], 4) ./ sqrt(sum(~isnan(lr_sum_ctrl), 4));
sem_lr_diff_exp  = nanstd(abs(lr_diff_exp),  [], 4) ./ sqrt(sum(~isnan(lr_diff_exp),  4));
sem_lr_sum_exp   = nanstd(abs(lr_sum_exp),   [], 4) ./ sqrt(sum(~isnan(lr_sum_exp),   4));

% Avoid division by zero
sem_lr_diff_ctrl(sem_lr_diff_ctrl==0) = NaN;
sem_lr_sum_ctrl(sem_lr_sum_ctrl==0)   = NaN;
sem_lr_diff_exp(sem_lr_diff_exp==0)   = NaN;
sem_lr_sum_exp(sem_lr_sum_exp==0)     = NaN;

% Within-group t-scores
t_lr_diff_ctrl = avg_lr_diff_ctrl ./ sem_lr_diff_ctrl;
t_lr_sum_ctrl  = avg_lr_sum_ctrl  ./ sem_lr_sum_ctrl;
t_lr_diff_exp  = avg_lr_diff_exp  ./ sem_lr_diff_exp;
t_lr_sum_exp   = avg_lr_sum_exp   ./ sem_lr_sum_exp;

% Pooled SEM for between-group difference
sem_diff_lr_diff = sqrt(sem_lr_diff_ctrl.^2 + sem_lr_diff_exp.^2);
sem_diff_lr_sum  = sqrt(sem_lr_sum_ctrl.^2  + sem_lr_sum_exp.^2);
sem_diff_lr_diff(isnan(sem_diff_lr_diff) | sem_diff_lr_diff==0) = NaN;
sem_diff_lr_sum(isnan(sem_diff_lr_sum)  | sem_diff_lr_sum==0)  = NaN;

% Between-group t-scores (Welch)
t_lr_diff_groupdiff = avg_lr_diff_groupdiff ./ sem_diff_lr_diff;
t_lr_sum_groupdiff  = avg_lr_sum_groupdiff  ./ sem_diff_lr_sum;

% Generate t-score groupdiff video
t_lim = [-6 6];
write_lr_video(t_lr_diff_groupdiff, t_lr_sum_groupdiff, half_atlas, brainMask_cropped, ...
    base_root, ['t_lr_diff_sum_groupdiff_' ctrl_type '_' exp_type '.mp4'], t_lim, ...
    [exp_type ' - ' ctrl_type ' (t-score)'], ...
    'LR abs difference t-score groupdiff', 'LR abs sum t-score groupdiff');

%% Plot surprise and surprise-masked t-scored videos

% Get number of mice per group
n_ctrl = sum(~isnan(lr_diff_ctrl), 4);
n_exp  = sum(~isnan(lr_diff_exp),  4);

% Compute Welch–Satterthwaite degrees of freedom
var1_diff = sem_lr_diff_ctrl.^2;   var2_diff = sem_lr_diff_exp.^2;
var1_sum  = sem_lr_sum_ctrl.^2;    var2_sum  = sem_lr_sum_exp.^2;
df_diff = (var1_diff + var2_diff).^2 ./ (var1_diff.^2./(n_ctrl-1) + var2_diff.^2./(n_exp-1));
df_sum  = (var1_sum  + var2_sum ).^2 ./ (var1_sum.^2 ./ (n_ctrl-1) + var2_sum.^2 ./ (n_exp-1));
df_diff(n_ctrl < 2 | n_exp < 2) = NaN;
df_sum(n_ctrl  < 2 | n_exp < 2) = NaN;

% Compute t-zest p-values and surprise
p_diff = 2 * tcdf(-abs(t_lr_diff_groupdiff), df_diff);
p_sum  = 2 * tcdf(-abs(t_lr_sum_groupdiff),  df_sum);
surp_diff = -log10(p_diff);
surp_sum  = -log10(p_sum);

% Generate surprise video
write_lr_video(surp_diff, surp_sum, half_atlas, brainMask_cropped, ...
    base_root, ['surp_lr_diff_sum_groupdiff_' ctrl_type '_' exp_type '.mp4'], [0 8], ...
    ['-log_{10}(p) | ' exp_type ' vs ' ctrl_type], ...
    'LR abs difference surprise', 'LR abs sum surprise');

% Generate surprise-masked t-score video
surp_thresh = -log10(0.05);
write_lr_video_surpmask( ...
    t_lr_diff_groupdiff, t_lr_sum_groupdiff, ...
    half_atlas, brainMask_cropped, ...
    base_root, ['t_lr_diff_sum_groupdiff_' ctrl_type '_' exp_type '_surpmask.mp4'], ...
    t_lim, [exp_type ' - ' ctrl_type ' (t-score, p<0.05)'], ...
    'LR abs diff t-score (masked)', 'LR abs sum t-score (masked)', ...
    surp_diff, surp_thresh);

% Display output
disp('Generalized LR analysis completed successfully!');