function plot_rescaling_test(run_settings)
%PLOT_RESCALING_TEST  The cluster test at each rescaling of the experimental mice.
%   PLOT_RESCALING_TEST(run_settings) does the work of run_rescaling_figure,
%   which sets the fields of run_settings and says what each one does.
%
%   Reads, from comp_out_dir, what run_mouse_influence and run_scale_free_test
%   wrote for the same comparison and smoothing, <tag> being
%   <comp_tag>_smooth<sigma>:
%     Influence_Slopes_<tag>.csv     all the mice at each slope of the
%                                    alignment, the fitted one included: each
%                                    sign's heaviest cluster in the region, its
%                                    voxels and one-sided p, and step 3's p of
%                                    either sign with the heavier sign
%     Influence_Folds_<tag>.csv      the slope with all the mice (fold 0) and
%                                    refitted without each mouse
%     Scale_Free_Clusters_<tag>.csv  the region's heaviest cluster of each sign
%                                    on the asymmetry index, and step 3's
%   Writes into comp_out_dir Rescaling_Test_<tag>, .fig and .png.
%
%   The slope is that of the line step 3 fits from the experimental group's
%   mean tissue profile along AP onto the control group's. On the test's maps
%   it multiplies every experimental mouse's |L - R| (the intercept cancels in
%   L - R, the common factor in the t); at 1 each group keeps its own scale of
%   step 2. Two panels against it. The p over every split of the mice: step
%   3's p, that of the heavier sign against each split's heavier (the p
%   quoted), and each sign's own, one-sided (for the experimental-higher sign
%   the test the scale-free test named in advance). And each sign's cluster
%   voxels. The lines join the slopes of the sweep; the fitted slope's values
%   are larger markers, edged in black, and step 3's p is joined only between
%   slopes with the same heavier sign. Both panels mark the fitted slope with
%   its jackknife standard error, from the slopes s_i refitted without each of
%   the n mice,
%       se = sqrt((n - 1) / n * sum((s_i - mean(s)).^2))
%   and the slope refitted without the mouse that moves it most; beside each
%   axis, the same p and clusters on the index |L - R| / (L + R) of
%   run_scale_free_test, where no scale enters. Under the title, the slopes
%   between which step 3's p crosses 0.05 for an experimental-higher cluster
%   and the heavier sign changes. Every number is read from the tables. The
%   fitted slope must be fold 0's, and its cluster step 3's of the scale-free
%   test's table, so the two steps read the same step 3.

% settings of run_rescaling_figure, under the names the code below uses
ctrl_type = run_settings.ctrl_type;
exp_type = run_settings.exp_type;
cluster_region = run_settings.cluster_region;
apply_smoothing = run_settings.apply_smoothing;
smooth_sigma = run_settings.smooth_sigma;
comp_tag = run_settings.comp_tag;
comp_out_dir = run_settings.comp_out_dir;

% the file names carry the smoothing, as run_mouse_influence's
if apply_smoothing
    smooth_suffix = sprintf('_smooth%g', smooth_sigma);
else
    smooth_suffix = '_nosmooth';
end
file_tag = [comp_tag smooth_suffix];

%% Read the tables

% all the mice at each slope, in order of the slope; the folds; the clusters of
% the scale-free test
slopes_file = fullfile(comp_out_dir, ['Influence_Slopes_' file_tag '.csv']);
folds_file = fullfile(comp_out_dir, ['Influence_Folds_' file_tag '.csv']);
free_file = fullfile(comp_out_dir, ['Scale_Free_Clusters_' file_tag '.csv']);
T_slopes = read_step_table(slopes_file, 'run_mouse_influence');
T_slopes = sortrows(T_slopes, 'slope');
T_folds = read_step_table(folds_file, 'run_mouse_influence');
T_free = read_step_table(free_file, 'run_scale_free_test');

% the fitted slope must be fold 0's, and its cluster step 3's of the scale-free
% test's table
fitted = T_slopes(T_slopes.is_fitted == 1, :);
fold_zero = T_folds(T_folds.fold == 0, :);
step3 = cluster_row(T_free, [exp_type ' higher, step 3 (test maps)'], free_file);
check_same_step3(fitted, fold_zero, step3, slopes_file, folds_file, free_file);

%% The slope's marks and the index's clusters

% the slope refitted without each mouse, its jackknife standard error, and the
% mouse whose leaving out moves it most
marks = slope_marks(T_folds, fitted.slope);

% the smallest p possible: one split for one sign; for either sign, two with
% groups of equal size, a split and its mirror being as heavy
floors = struct();
floors.one_sign = 1 / fitted.n_splits;
floors.either_sign = (1 + (fold_zero.n_ctrl == fold_zero.n_exp)) / fitted.n_splits;

% the region's two clusters on the index, where no scale enters
free = struct();
free.exp = cluster_row(T_free, [exp_type ' higher, index'], free_file);
free.ctrl = cluster_row(T_free, [ctrl_type ' higher, index'], free_file);

%% Figure

comparison = struct('ctrl_type', ctrl_type, 'exp_type', exp_type, ...
    'cluster_region', cluster_region);
print_summary(T_slopes, marks, free, comparison);
plot_figure(T_slopes, marks, floors, free, comparison, file_tag, comp_out_dir);
fprintf('Rescaling figure saved to: %s\n', ...
    fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.png']));

end

% ===== Local functions: reading =====

function T = read_step_table(file, step)
% A table written by an earlier step, its text columns as char.

if ~isfile(file)
    error(['plot_rescaling_test: %s not found. It is written by %s with the same ' ...
           'comparison and smoothing; run that first.'], file, step);
end
T = readtable(file, 'TextType', 'char');
end

function row = cluster_row(T_free, name, free_file)
% The row of the scale-free test's clusters named name.

row = T_free(strcmp(T_free.cluster, name), :);
if height(row) ~= 1
    error('plot_rescaling_test: no single row ''%s'' in %s (its rows: %s).', name, ...
        free_file, strjoin(T_free.cluster', '; '));
end
end

function check_same_step3(fitted, fold_zero, step3, slopes_file, folds_file, free_file)
% The fitted slope of the slopes' table is fold 0's, and its cluster is step
% 3's of the scale-free test's table: voxels, mass, one-sided p and the p of
% either sign.

% the tables hold 15 significant digits
tolerance = 1e-9;
if height(fitted) ~= 1 || height(fold_zero) ~= 1
    error(['plot_rescaling_test: expected one fitted slope in %s and one fold 0 ' ...
           'in %s, found %d and %d.'], slopes_file, folds_file, height(fitted), ...
           height(fold_zero));
end
if abs(fitted.slope - fold_zero.slope) > tolerance
    error(['plot_rescaling_test: the fitted slope of %s (%.6f) is not fold 0''s ' ...
           'of %s (%.6f). Rerun run_mouse_influence.'], slopes_file, fitted.slope, ...
           folds_file, fold_zero.slope);
end
is_same = fitted.cluster_n == step3.voxels && ...
    abs(fitted.cluster_mass - step3.mass) <= tolerance * step3.mass && ...
    abs(fitted.p_positive - step3.p_one_sided) <= tolerance && ...
    abs(fitted.p_either_sign - step3.p_either_sign) <= tolerance;
if ~is_same
    error(['plot_rescaling_test: the cluster at the fitted slope in %s (%d voxels, ' ...
           'mass %.2f, p %.4f, either sign %.4f) is not step 3''s cluster in %s ' ...
           '(%d voxels, mass %.2f, p %.4f, either sign %.4f): the two steps did not ' ...
           'read the same step 3. Rerun run_mouse_influence and run_scale_free_test ' ...
           'on the same data.'], slopes_file, fitted.cluster_n, fitted.cluster_mass, ...
           fitted.p_positive, fitted.p_either_sign, free_file, step3.voxels, ...
           step3.mass, step3.p_one_sided, step3.p_either_sign);
end
end

function marks = slope_marks(T_folds, fitted_slope)
% The fitted slope; the slopes refitted without each mouse and their jackknife
% standard error; the mouse whose leaving out moves the slope most, its slope
% and its group.

is_fold = T_folds.fold > 0;
refitted = T_folds.slope(is_fold);
left_out = T_folds.left_out(is_fold);
groups = T_folds.left_out_group(is_fold);
n_mice = numel(refitted);
[~, i_most] = max(abs(refitted - fitted_slope));

marks = struct();
marks.fitted = fitted_slope;
marks.se = sqrt((n_mice - 1) / n_mice * sum((refitted - mean(refitted)).^2));
marks.n_mice = n_mice;
marks.without = refitted(i_most);
marks.without_mouse = strrep(left_out{i_most}, '_Gria1', '');
marks.without_group = groups{i_most};
end

function lines = takeaway_line(T_slopes, comparison)
% The two lines under the title, from the slopes' table: step 3's p at the
% fitted slope; where it crosses 0.05 for an experimental-higher cluster, and
% where the heavier sign is the control group's.

fitted = T_slopes(T_slopes.is_fitted == 1, :);
is_exp_below = T_slopes.either_sign > 0 & T_slopes.p_either_sign < 0.05;
is_ctrl_heavier = T_slopes.either_sign < 0;
lines = {
    sprintf(['Step 3''s p, %.3f at the fitted slope %.2f, depends on the rescaling ' ...
             'of the %s mice:'], fitted.p_either_sign, fitted.slope, comparison.exp_type)
    sprintf(['for the %s-higher cluster it is under 0.05 %s, and the %s-higher ' ...
             'cluster is the heavier %s.'], comparison.exp_type, ...
             slopes_where(T_slopes.slope, is_exp_below), comparison.ctrl_type, ...
             slopes_where(T_slopes.slope, is_ctrl_heavier))
    };
end

function words = slopes_where(slopes, is_true)
% The slopes where is_true holds, in words: from or up to a slope between two
% of them when it changes once, else the slopes themselves.

i_change = find(diff(is_true) ~= 0);
if ~any(is_true)
    words = 'at no slope';
elseif all(is_true)
    words = 'at every slope';
elseif isscalar(i_change) && is_true(end)
    words = sprintf('from a slope between %.2f and %.2f on', slopes(i_change), ...
        slopes(i_change + 1));
elseif isscalar(i_change)
    words = sprintf('up to a slope between %.2f and %.2f', slopes(i_change), ...
        slopes(i_change + 1));
else
    words = ['at slopes ' strjoin(compose('%.2f', slopes(is_true)), ', ')];
end
end

function print_summary(T_slopes, marks, free, comparison)
% The numbers the figure shows.

fprintf('%s, %s against %s, the heaviest cluster of each sign over %d splits:\n', ...
    comparison.cluster_region, comparison.exp_type, comparison.ctrl_type, ...
    T_slopes.n_splits(1));
for s = 1:height(T_slopes)
    if T_slopes.is_fitted(s) == 1
        note = ' (fitted)';
    else
        note = '';
    end
    fprintf(['  slope %.2f%s: %s higher %d voxels, p %.3f; %s higher %d voxels, ' ...
             'p %.3f; either sign p %.3f, %s higher the heavier\n'], ...
             T_slopes.slope(s), note, comparison.exp_type, T_slopes.cluster_n(s), ...
             T_slopes.p_positive(s), comparison.ctrl_type, T_slopes.negative_n(s), ...
             T_slopes.p_negative(s), T_slopes.p_either_sign(s), ...
             heavier_name(T_slopes.either_sign(s), comparison));
end
fprintf(['  fitted slope %.4f, jackknife standard error %.3f over %d mice; refitted ' ...
         'without %s %.4f\n'], marks.fitted, marks.se, marks.n_mice, ...
         marks.without_mouse, marks.without);
fprintf(['  index: %s higher %d voxels, p %.3f; %s higher %d voxels, p %.3f; ' ...
         'either sign p %.3f\n'], comparison.exp_type, free.exp.voxels, ...
         free.exp.p_one_sided, comparison.ctrl_type, free.ctrl.voxels, ...
         free.ctrl.p_one_sided, free_either_sign(free));
takeaway = takeaway_line(T_slopes, comparison);
fprintf('  %s %s\n', takeaway{:});
end

function name = heavier_name(either_sign, comparison)
% The group whose higher cluster is the heavier, from the sign of step 3's score.

name = comparison.exp_type;
if either_sign < 0
    name = comparison.ctrl_type;
end
end

function [p, either_sign] = free_either_sign(free)
% The index's p of either sign, given on the row of the heavier sign only.

either_sign = 1;
p = free.exp.p_either_sign;
if isnan(p)
    either_sign = -1;
    p = free.ctrl.p_either_sign;
end
end

% ===== Local functions: figure =====

function plot_figure(T_slopes, marks, floors, free, comparison, file_tag, comp_out_dir)
% The p and the size of each sign's cluster against the slope, each panel with
% the index's values on a narrow axis beside it; the title, the line to take
% from it, the legend and a note on the readings.

% the series' colours: each sign's as in run_mouse_influence's slope panels,
% step 3's p of either sign in black, the one quoted
colours = struct();
colours.exp = sep_palette('experimental_mean');
colours.ctrl = sep_palette('control_mean');
colours.either = [0 0 0];
x_limits = [min(T_slopes.slope) - 0.05, max(T_slopes.slope) + 0.05];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

% the p, and the index's beside it
axes('Position', [0.06 0.24 0.31 0.58]);
p_limits = [2e-3 1.2];
p_ticks = [0.005 0.01 0.02 0.05 0.1 0.2 0.5 1];
handles = draw_p_panel(T_slopes, marks, floors, colours, x_limits, p_limits, ...
    p_ticks, comparison);
axes('Position', [0.38 0.24 0.05 0.58]);
draw_free_side(free_p_points(free, colours), '%.3f', p_limits, p_ticks);

% the size, and the index's beside it
% (a margin under 0, so a sign without a cluster shows above the axis)
axes('Position', [0.535 0.24 0.31 0.58]);
size_top = 1.1 * max([T_slopes.cluster_n; T_slopes.negative_n; free.exp.voxels]);
size_limits = [-0.04 1] * size_top;
draw_size_panel(T_slopes, marks, colours, x_limits, size_limits, comparison);
axes('Position', [0.855 0.24 0.05 0.58]);
draw_free_side(free_size_points(free, colours), '%d', size_limits, []);

% the legend, under both panels
labels = {
    sprintf('%s higher, one-sided', comparison.exp_type)
    sprintf('%s higher, one-sided', comparison.ctrl_type)
    sprintf(['either sign: step 3''s p, the one quoted (square: the %s-higher ' ...
             'cluster the heavier)'], comparison.ctrl_type)
    sprintf('fitted slope %.2f (large markers)', marks.fitted)
    sprintf('\\pm 1 jackknife SE (%.2f, %d mice)', marks.se, marks.n_mice)
    sprintf('refitted without %s (%s), %.2f', marks.without_mouse, ...
        marks.without_group, marks.without)
    };
lgd = legend(handles, labels, 'NumColumns', 3, 'FontSize', 11, 'Box', 'off');
lgd.Units = 'normalized';
lgd.Position(1) = 0.5 - lgd.Position(3) / 2;
lgd.Position(2) = 0.085;

% a note on the readings
note = {
    sprintf(['Slopes: step 3''s maps, every %s mouse''s |L - R| times the slope; ' ...
             'the splits keep it (run_mouse_influence). No scale: step 3''s test on ' ...
             '|L - R| / (L + R) of the collected stacks, with neither step 2''s line ' ...
             'per mouse nor the alignment (run_scale_free_test).'], comparison.exp_type)
    sprintf(['One-sided: a sign''s heaviest cluster against the same sign''s in ' ...
             'every split (%s higher: the test the scale-free test named in ' ...
             'advance). Either sign: the heavier against each split''s heavier. ' ...
             'Dotted: the smallest p possible.'], comparison.exp_type)
    };
annotation('textbox', [0.04 0.01 0.92 0.05], 'String', note, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
    'Interpreter', 'none');

title_line = sprintf(['%s, %s against %s: the cluster test at each rescaling of the ' ...
                      '%s mice - %s'], comparison.cluster_region, comparison.exp_type, ...
                      comparison.ctrl_type, comparison.exp_type, ...
                      strrep(file_tag, '_', ' '));
takeaway = takeaway_line(T_slopes, comparison);
sgtitle({title_line, ['\rm\fontsize{13}' takeaway{1}], ...
    ['\rm\fontsize{13}' takeaway{2}]}, ...
    'FontSize', 15, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.png']), ...
    'Resolution', 300);
end

function handles = draw_p_panel(T_slopes, marks, floors, colours, x_limits, ...
    y_limits, y_ticks, comparison)
% Each sign's p and step 3's p at each slope on a log axis, the slope's marks
% behind them, a dashed line at 0.05 and dotted ones at the smallest p; the
% values at the fitted slope. Returns the handles the legend names.

hold on;
box on;
grid on;
mark_handles = draw_slope_marks(marks, y_limits, y_limits(2) / 1.4);
yline(0.05, '--', 'p = 0.05', 'Color', [0.4 0.4 0.4], 'FontSize', 10, ...
    'LabelHorizontalAlignment', 'left');
yline(floors.one_sign, ':', 'one split', 'Color', [0.5 0.5 0.5], 'FontSize', 10, ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');
yline(floors.either_sign, ':', 'either sign''s smallest', 'Color', [0.5 0.5 0.5], ...
    'FontSize', 10, 'LabelHorizontalAlignment', 'right', ...
    'LabelVerticalAlignment', 'bottom');

% each sign's p over the sweep's slopes, and at the fitted slope
sweep = T_slopes(T_slopes.is_fitted == 0, :);
fitted = T_slopes(T_slopes.is_fitted == 1, :);
exp_line = draw_series(sweep.slope, sweep.p_positive, 'o', colours.exp, colours.exp);
ctrl_line = draw_series(sweep.slope, sweep.p_negative, 'o', colours.ctrl, colours.ctrl);
draw_fitted_point(marks.fitted, fitted.p_positive, 'o', colours.exp);
draw_fitted_point(marks.fitted, fitted.p_negative, 'o', colours.ctrl);

% step 3's p, joined only between slopes with the same heavier sign: across the
% change the p is another cluster's; a square where the control-higher one is
% the heavier
run_id = cumsum([1; diff(sweep.either_sign) ~= 0]);
for r = 1:run_id(end)
    in_run = run_id == r;
    marker = sign_marker(sweep.either_sign(find(in_run, 1)));
    either_line = draw_series(sweep.slope(in_run), sweep.p_either_sign(in_run), ...
        marker, colours.either, 'w');
end
draw_fitted_point(marks.fitted, fitted.p_either_sign, ...
    sign_marker(fitted.either_sign), 'w');

% the values at the fitted slope, left of its line, above the curves
value_lines = {
    sprintf('at the fitted slope, %.2f:', marks.fitted), [0.3 0.3 0.3]
    sprintf('step 3''s p, either sign %.3f', fitted.p_either_sign), colours.either
    sprintf('%s higher, one-sided %.3f', comparison.exp_type, fitted.p_positive), ...
        colours.exp
    sprintf('%s higher, one-sided %.3f', comparison.ctrl_type, fitted.p_negative), ...
        colours.ctrl
    };
for k = 1:size(value_lines, 1)
    text(marks.fitted - 0.03, 0.5 / 1.35^(k - 1), value_lines{k, 1}, ...
        'HorizontalAlignment', 'right', 'Color', value_lines{k, 2}, 'FontSize', 10);
end

set(gca, 'YScale', 'log', 'FontSize', 11, 'Layer', 'top');
ylim(y_limits);
yticks(y_ticks);
xlim(x_limits);
xlabel(sprintf('rescaling of the %s mice: slope of step 3''s alignment (1: none)', ...
    comparison.exp_type), 'FontSize', 11);
ylabel(sprintf('p over every split of the mice (%d)', fitted.n_splits), 'FontSize', 11);
title({sprintf('p of the cluster mass in %s', comparison.cluster_region), ...
    ['\rm\fontsize{10}the heaviest cluster where |L - R| is higher in each group, ' ...
    'and step 3''s p'], ''}, 'FontSize', 13);
handles = [exp_line, ctrl_line, either_line, mark_handles];
end

function draw_size_panel(T_slopes, marks, colours, x_limits, y_limits, comparison)
% Each sign's cluster voxels at each slope, the slope's marks behind them, and
% the values at the fitted slope.

hold on;
box on;
grid on;
draw_slope_marks(marks, y_limits, y_limits(2) - 0.04 * diff(y_limits));
sweep = T_slopes(T_slopes.is_fitted == 0, :);
fitted = T_slopes(T_slopes.is_fitted == 1, :);
draw_series(sweep.slope, sweep.cluster_n, 'o', colours.exp, colours.exp);
draw_series(sweep.slope, sweep.negative_n, 'o', colours.ctrl, colours.ctrl);
draw_fitted_point(marks.fitted, fitted.cluster_n, 'o', colours.exp);
draw_fitted_point(marks.fitted, fitted.negative_n, 'o', colours.ctrl);

% the values at the fitted slope, left of its line, above the curves
value_lines = {
    sprintf('at the fitted slope, %.2f:', marks.fitted), [0.3 0.3 0.3]
    sprintf('%s higher %d voxels', comparison.exp_type, fitted.cluster_n), colours.exp
    sprintf('%s higher %d voxels', comparison.ctrl_type, fitted.negative_n), ...
        colours.ctrl
    };
for k = 1:size(value_lines, 1)
    text(marks.fitted - 0.03, (0.86 - 0.05 * (k - 1)) * y_limits(2), ...
        value_lines{k, 1}, 'HorizontalAlignment', 'right', 'Color', ...
        value_lines{k, 2}, 'FontSize', 10);
end

set(gca, 'FontSize', 11, 'Layer', 'top');
ylim(y_limits);
xlim(x_limits);
xlabel(sprintf('rescaling of the %s mice: slope of step 3''s alignment (1: none)', ...
    comparison.exp_type), 'FontSize', 11);
ylabel('voxels of 10 um (1,000 voxels = 0.001 mm^3)', 'FontSize', 11);
title({sprintf('Size of each cluster in %s', comparison.cluster_region), ...
    '\rm\fontsize{10}the same clusters', ''}, 'FontSize', 13);
end

function handle = draw_series(x, y, marker, colour, face)
% One series joined by a line, its markers in face.

handle = plot(x, y, ['-' marker], 'Color', colour, 'MarkerFaceColor', face, ...
    'LineWidth', 1.6, 'MarkerSize', 6);
end

function draw_fitted_point(x, y, marker, face)
% A value at the fitted slope: a larger marker edged in black.

plot(x, y, marker, 'MarkerFaceColor', face, 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.4, 'MarkerSize', 10);
end

function marker = sign_marker(either_sign)
% Step 3's p: a circle where the experimental-higher cluster is the heavier, a
% square where the control-higher one is.

marker = 'o';
if either_sign < 0
    marker = 's';
end
end

function handles = draw_slope_marks(marks, y_limits, label_y)
% The fitted slope with a band of +/- 1 jackknife standard error, and the slope
% refitted without the mouse that moves it most, over the panel's height, each
% line named at label_y. Returns the fitted line, the band and the refitted line.

band = patch(marks.fitted + marks.se * [-1 1 1 -1], y_limits([1 1 2 2]), ...
    sep_palette('paired_lines'), 'FaceAlpha', 0.15, 'EdgeColor', ...
    sep_palette('paired_lines'), 'LineWidth', 0.5);
fitted_line = plot([1 1] * marks.fitted, y_limits, '-', 'Color', [0.3 0.3 0.3], ...
    'LineWidth', 1.2);
without_line = plot([1 1] * marks.without, y_limits, ':', 'Color', [0.3 0.3 0.3], ...
    'LineWidth', 2);
text(marks.fitted + 0.01, label_y, sprintf('fitted %.2f \\pm %.2f', marks.fitted, ...
    marks.se), 'Color', [0.3 0.3 0.3], 'FontSize', 10);
text(marks.without + 0.01, label_y, sprintf('without %s %.2f', marks.without_mouse, ...
    marks.without), 'Color', [0.3 0.3 0.3], 'FontSize', 10);
handles = [fitted_line, band, without_line];
end

function points = free_p_points(free, colours)
% The index's p beside the p panel, in the panel's markers: each sign's
% one-sided p, and the p of either sign.

[p_either, either_sign] = free_either_sign(free);
points = struct( ...
    'value', {free.exp.p_one_sided, free.ctrl.p_one_sided, p_either}, ...
    'is_none', {free.exp.voxels == 0, free.ctrl.voxels == 0, false}, ...
    'colour', {colours.exp, colours.ctrl, colours.either}, ...
    'face', {colours.exp, colours.ctrl, [1 1 1]}, ...
    'marker', {'o', 'o', sign_marker(either_sign)});
end

function points = free_size_points(free, colours)
% The index's cluster voxels beside the size panel, each sign's.

points = struct( ...
    'value', {free.exp.voxels, free.ctrl.voxels}, ...
    'is_none', {free.exp.voxels == 0, free.ctrl.voxels == 0}, ...
    'colour', {colours.exp, colours.ctrl}, ...
    'face', {colours.exp, colours.ctrl}, ...
    'marker', {'o', 'o'});
end

function draw_free_side(points, value_format, y_limits, y_ticks)
% The region's values on the index, beside a panel's axis and on its scale,
% each with its value, or open and 'none' where the sign has no cluster (its
% p 1, its voxels 0); a log axis with 0.05 dashed when y_ticks are given.

hold on;
box on;
grid on;
for k = 1:numel(points)
    face = points(k).face;
    value_text = sprintf(value_format, points(k).value);
    if points(k).is_none
        face = [1 1 1];
        value_text = 'none';
    end
    plot(0.3, points(k).value, points(k).marker, 'Color', points(k).colour, ...
        'MarkerFaceColor', face, 'LineWidth', 1.4, 'MarkerSize', 8);
    text(0.45, points(k).value, value_text, 'Color', points(k).colour, 'FontSize', 10);
end
if ~isempty(y_ticks)
    set(gca, 'YScale', 'log');
    yticks(y_ticks);
    yline(0.05, '--', 'Color', [0.4 0.4 0.4]);
end
set(gca, 'FontSize', 11, 'YTickLabel', [], 'XTick', []);
ylim(y_limits);
xlim([0 1]);
xlabel({'no scale', '|L - R| / (L + R)'}, 'FontSize', 11);
end
