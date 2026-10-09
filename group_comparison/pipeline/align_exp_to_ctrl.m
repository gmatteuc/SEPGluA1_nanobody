function [interest_region, norm_ctrl, norm_exp, slope, intercept, norm_ctrl_med_fact, ...
    norm_exp_med_fact] = align_exp_to_ctrl(med_data_4d_ctrl, med_data_4d_exp)
%ALIGN_EXP_TO_CTRL  The line from the experimental group's profile to the control's.
%   [interest_region, norm_ctrl, norm_exp, slope, intercept, norm_ctrl_med_fact,
%   norm_exp_med_fact] = ALIGN_EXP_TO_CTRL(med_data_4d_ctrl, med_data_4d_exp)
%   fits a line from the experimental group's mean plane profile onto the
%   control group's over planes 200 to 700, from the profiles of
%   plane_tissue_means (planes x mice each). Returns the fit planes, both
%   groups' profiles after the alignment, the line, and the common scale
%   factor of both groups (the two outputs are the same number). A volume of
%   the experimental group goes onto the common scale as
%       ((volume .* slope) + intercept) ./ factor
%   one of the control group as volume ./ factor. Used by group_differences
%   and per_mouse_region_values.

% planes the line is fitted on: 200 to 700 of the 900 (379 to 879 of the 10 um
% annotation), away from both ends of the brain
interest_region = 200:700;

% the groups' mean plane profiles
mean_profile_ctrl = nanmean(med_data_4d_ctrl, 2); %#ok<*NANMEAN>
mean_profile_exp = nanmean(med_data_4d_exp, 2);

% the line from the experimental profile to the control one, on the fit planes
% where both have a value
y_target = mean_profile_ctrl(interest_region);
x_source = mean_profile_exp(interest_region);
valid_idx = ~isnan(x_source) & ~isnan(y_target);
x_source = x_source(valid_idx);
y_target = y_target(valid_idx);
p = polyfit(x_source, y_target, 1);
slope = p(1);
intercept = p(2);

fprintf('Alignment Parameters (Exp -> Ctrl): Slope = %.4f, Intercept = %.4f\n', slope, ...
    intercept);

% the control profiles stay as they are, the experimental ones go through the line
norm_ctrl = med_data_4d_ctrl;
norm_exp = (med_data_4d_exp .* slope) + intercept;

% one common scale for both groups: the average of the two groups' mean
% intensity over planes 300-500, after the alignment, so every map below is in
% units of the mid-brain tissue mean
interest_region_bis = 300:500;

% each group's mean over the planes, then over its mice
norm_ctrl_med_fact = nanmean(nanmean(med_data_4d_ctrl(interest_region_bis, :), 1));
norm_exp_med_fact = nanmean(nanmean(norm_exp(interest_region_bis, :), 1));

% one factor for both: it sets the unit of the maps and leaves every t unchanged
avg_med_fact = (norm_ctrl_med_fact + norm_exp_med_fact)./2;
norm_ctrl_med_fact = avg_med_fact;
norm_exp_med_fact = avg_med_fact;
end
