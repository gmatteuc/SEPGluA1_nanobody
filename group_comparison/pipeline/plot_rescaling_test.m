function plot_rescaling_test(run_settings)
%PLOT_RESCALING_TEST  The region's cluster test at each rescaling of the experimental mice.
%   PLOT_RESCALING_TEST(run_settings) does the work of run_rescaling_figure,
%   which sets the fields of run_settings and says what each one does.
%
%   Reads, from comp_out_dir, what run_mouse_influence and run_scale_free_test
%   wrote for the same comparison and smoothing, <tag> being
%   <comp_tag>_smooth<sigma>:
%     Influence_Slopes_<tag>.csv     all the mice at each slope of the
%                                    alignment, the fitted one included: each
%                                    sign's heaviest cluster in the region, its
%                                    voxels and p, and which sign is the heavier
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
%   step 2. Two panels against it: the p of each sign's heaviest cluster mass
%   over every split of the mice (for the experimental-higher sign one-sided,
%   the test named in advance), and each cluster's voxels. Both mark the
%   fitted slope with its jackknife standard error, from the slopes s_i
%   refitted without each of the n mice,
%       se = sqrt((n - 1) / n * sum((s_i - mean(s)).^2))
%   and the slope refitted without the mouse that moves it most; beside each
%   axis, the same two clusters on the index |L - R| / (L + R) of
%   run_scale_free_test, where no scale enters. Under the title, the slopes
%   between which the heavier sign changes. Every number is read from the
%   tables. The fitted slope must be fold 0's, and its cluster step 3's of the
%   scale-free test's table, so the two steps read the same step 3.

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

% the region's two clusters on the index, where no scale enters
free = struct();
free.exp = cluster_row(T_free, [exp_type ' higher, index'], free_file);
free.ctrl = cluster_row(T_free, [ctrl_type ' higher, index'], free_file);

%% Figure

comparison = struct('ctrl_type', ctrl_type, 'exp_type', exp_type, ...
    'cluster_region', cluster_region);
print_summary(T_slopes, marks, free, comparison);
plot_figure(T_slopes, marks, free, comparison, file_tag, comp_out_dir);
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
% 3's of the scale-free test's table: voxels, mass and one-sided p.

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
    abs(fitted.p_positive - step3.p_one_sided) <= tolerance;
if ~is_same
    error(['plot_rescaling_test: the cluster at the fitted slope in %s (%d voxels, ' ...
           'mass %.2f, p %.4f) is not step 3''s cluster in %s (%d voxels, mass %.2f, ' ...
           'p %.4f): the two steps did not read the same step 3. Rerun ' ...
           'run_mouse_influence and run_scale_free_test on the same data.'], ...
           slopes_file, fitted.cluster_n, fitted.cluster_mass, fitted.p_positive, ...
           free_file, step3.voxels, step3.mass, step3.p_one_sided);
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

function line = heavier_line(T_slopes, comparison)
% The line under the title: the slopes between which the heavier of the two
% clusters changes sign, from the slopes' table.

ctrl_heavier = T_slopes.either_sign < 0;
i_last = find(ctrl_heavier, 1, 'last');
start = ['The p of the test named in advance depends on the rescaling between the ' ...
         'groups: '];
if isempty(i_last)
    line = sprintf('%sthe %s-higher patch is the heavier at every slope.', start, ...
        comparison.exp_type);
elseif all(ctrl_heavier(1:i_last)) && i_last < height(T_slopes)
    line = sprintf(['%sbelow a slope between %g and %g the %s-higher patch is the ' ...
                    'heavier.'], start, T_slopes.slope(i_last), ...
                    T_slopes.slope(i_last + 1), comparison.ctrl_type);
else
    slopes_text = strjoin(arrayfun(@(s) sprintf('%.2f', s), ...
        T_slopes.slope(ctrl_heavier)', 'UniformOutput', false), ', ');
    line = sprintf('%sthe %s-higher patch is the heavier at slopes %s.', start, ...
        comparison.ctrl_type, slopes_text);
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
             'p %.3f; either sign p %.3f\n'], T_slopes.slope(s), note, ...
             comparison.exp_type, T_slopes.cluster_n(s), T_slopes.p_positive(s), ...
             comparison.ctrl_type, T_slopes.negative_n(s), T_slopes.p_negative(s), ...
             T_slopes.p_either_sign(s));
end
fprintf(['  fitted slope %.4f, jackknife standard error %.3f over %d mice; refitted ' ...
         'without %s %.4f\n'], marks.fitted, marks.se, marks.n_mice, ...
         marks.without_mouse, marks.without);
fprintf('  index: %s higher %d voxels, p %.3f; %s higher %d voxels, p %.3f\n', ...
    comparison.exp_type, free.exp.voxels, free.exp.p_one_sided, comparison.ctrl_type, ...
    free.ctrl.voxels, free.ctrl.p_one_sided);
fprintf('  %s\n', heavier_line(T_slopes, comparison));
end

% ===== Local functions: figure =====

function plot_figure(T_slopes, marks, free, comparison, file_tag, comp_out_dir)
% The p and the size of each sign's cluster against the slope, each panel with
% the index's clusters on a narrow axis beside it; the title, the line to take
% from it, the legend and a note on the two readings.

% the clusters' colours, as in run_mouse_influence's slope panels
colours = struct();
colours.exp = sep_palette('experimental_mean');
colours.ctrl = sep_palette('control_mean');
x_limits = [min(T_slopes.slope) - 0.05, max(T_slopes.slope) + 0.05];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

% the p, and the index's beside it
axes('Position', [0.06 0.21 0.31 0.61]);
p_limits = [2e-3 1.2];
handles = draw_p_panel(T_slopes, marks, colours, x_limits, p_limits, comparison);
axes('Position', [0.38 0.21 0.05 0.61]);
draw_free_side([free.exp.p_one_sided, free.ctrl.p_one_sided], ...
    [free.exp.voxels, free.ctrl.voxels], '%.3f', p_limits, true, colours);

% the size, and the index's beside it
% (a margin under 0, so a sign without a cluster shows above the axis)
axes('Position', [0.535 0.21 0.31 0.61]);
size_top = 1.1 * max([T_slopes.cluster_n; T_slopes.negative_n; free.exp.voxels]);
size_limits = [-0.04 1] * size_top;
draw_size_panel(T_slopes, marks, colours, x_limits, size_limits, comparison);
axes('Position', [0.855 0.21 0.05 0.61]);
draw_free_side([free.exp.voxels, free.ctrl.voxels], ...
    [free.exp.voxels, free.ctrl.voxels], '%d', size_limits, false, colours);

% the legend, under both panels
labels = {
    sprintf('%s higher (the test named in advance, one-sided)', comparison.exp_type)
    sprintf('%s higher', comparison.ctrl_type)
    sprintf('fitted slope, %.2f', marks.fitted)
    sprintf('\\pm 1 jackknife SE (%.2f, %d mice)', marks.se, marks.n_mice)
    sprintf('refitted without %s (%s), %.2f', marks.without_mouse, ...
        marks.without_group, marks.without)
    };
lgd = legend(handles, labels, 'Orientation', 'horizontal', 'NumColumns', 5, ...
    'FontSize', 11, 'Box', 'off');
lgd.Units = 'normalized';
lgd.Position(1) = 0.5 - lgd.Position(3) / 2;
lgd.Position(2) = 0.085;

% a note on the two readings
note = sprintf(['Slopes: step 3''s maps, every %s mouse''s |L - R| times the ' ...
                'slope; the splits keep it (run_mouse_influence). No scale: step ' ...
                '3''s test on |L - R| / (L + R) of the collected stacks, with ' ...
                'neither step 2''s line per mouse nor the alignment ' ...
                '(run_scale_free_test).'], comparison.exp_type);
annotation('textbox', [0.06 0.02 0.88 0.04], 'String', note, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
    'Interpreter', 'none');

title_line = sprintf(['%s, %s against %s: the cluster test at each rescaling of the ' ...
                      '%s mice - %s'], comparison.cluster_region, comparison.exp_type, ...
                      comparison.ctrl_type, comparison.exp_type, ...
                      strrep(file_tag, '_', ' '));
sgtitle({title_line, ['\rm\fontsize{13}' heavier_line(T_slopes, comparison)]}, ...
    'FontSize', 15, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.png']), ...
    'Resolution', 300);
end

function handles = draw_p_panel(T_slopes, marks, colours, x_limits, y_limits, comparison)
% Each sign's p at each slope on a log axis, the slope's marks behind them, a
% dashed line at 0.05 and a dotted one at one split; at the fitted slope its p
% and step 3's p of either sign. Returns the handles the legend names.

hold on;
box on;
grid on;
mark_handles = draw_slope_marks(marks, y_limits, y_limits(2) / 1.4);
exp_line = plot(T_slopes.slope, T_slopes.p_positive, '-o', 'Color', colours.exp, ...
    'MarkerFaceColor', colours.exp, 'LineWidth', 1.6, 'MarkerSize', 6);
ctrl_line = plot(T_slopes.slope, T_slopes.p_negative, '-o', 'Color', colours.ctrl, ...
    'MarkerFaceColor', colours.ctrl, 'LineWidth', 1.6, 'MarkerSize', 6);
yline(0.05, '--', 'p = 0.05', 'Color', [0.4 0.4 0.4], 'FontSize', 10, ...
    'LabelHorizontalAlignment', 'left');
n_splits = T_slopes.n_splits(1);
yline(1 / n_splits, ':', 'one split', 'Color', [0.5 0.5 0.5], 'FontSize', 10, ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');

% the fitted slope's p, below and to the left of its point, where no line runs
fitted = T_slopes.is_fitted == 1;
text(marks.fitted - 0.02, T_slopes.p_positive(fitted) / 1.4, ...
    sprintf('fitted: p %.3f (either sign, step 3''s: %.3f)', ...
    T_slopes.p_positive(fitted), T_slopes.p_either_sign(fitted)), ...
    'HorizontalAlignment', 'right', ...
    'Color', colours.exp, 'FontSize', 10);

set(gca, 'YScale', 'log', 'FontSize', 11, 'Layer', 'top');
ylim(y_limits);
yticks([0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]);
xlim(x_limits);
xlabel(sprintf('rescaling of the %s mice: slope of step 3''s alignment (1: none)', ...
    comparison.exp_type), 'FontSize', 11);
ylabel(sprintf('p over every split of the mice (%d)', n_splits), 'FontSize', 11);
title({sprintf('p of each cluster''s mass in %s', comparison.cluster_region), ...
    sprintf(['\\rm\\fontsize{10}the heaviest cluster where |L - R| is higher in each ' ...
    'group, one-sided']), ''}, 'FontSize', 13);
handles = [exp_line, ctrl_line, mark_handles];
end

function draw_size_panel(T_slopes, marks, colours, x_limits, y_limits, comparison)
% Each sign's cluster voxels at each slope, the slope's marks behind them, and
% the fitted slope's cluster size.

hold on;
box on;
grid on;
draw_slope_marks(marks, y_limits, y_limits(2) - 0.04 * diff(y_limits));
plot(T_slopes.slope, T_slopes.cluster_n, '-o', 'Color', colours.exp, ...
    'MarkerFaceColor', colours.exp, 'LineWidth', 1.6, 'MarkerSize', 6);
plot(T_slopes.slope, T_slopes.negative_n, '-o', 'Color', colours.ctrl, ...
    'MarkerFaceColor', colours.ctrl, 'LineWidth', 1.6, 'MarkerSize', 6);

% the fitted slope's cluster size, above and to the left of its point
fitted = T_slopes.is_fitted == 1;
text(marks.fitted - 0.02, T_slopes.cluster_n(fitted) + 0.05 * y_limits(2), ...
    sprintf('fitted: %d voxels', T_slopes.cluster_n(fitted)), ...
    'HorizontalAlignment', 'right', 'Color', colours.exp, 'FontSize', 10);

set(gca, 'FontSize', 11, 'Layer', 'top');
ylim(y_limits);
xlim(x_limits);
xlabel(sprintf('rescaling of the %s mice: slope of step 3''s alignment (1: none)', ...
    comparison.exp_type), 'FontSize', 11);
ylabel('voxels of 10 um (1,000 voxels = 0.001 mm^3)', 'FontSize', 11);
title({sprintf('Size of each cluster in %s', comparison.cluster_region), ...
    '\rm\fontsize{10}the same clusters', ''}, 'FontSize', 13);
end

function handles = draw_slope_marks(marks, y_limits, label_y)
% The fitted slope with a band of +/- 1 jackknife standard error, and the slope
% refitted without the mouse that moves it most, over the panel's height, each
% line named at label_y. Returns the fitted line, the band and the refitted line.

band = patch(marks.fitted + marks.se * [-1 1 1 -1], y_limits([1 1 2 2]), ...
    sep_palette('paired_lines'), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
fitted_line = plot([1 1] * marks.fitted, y_limits, '-', 'Color', [0.3 0.3 0.3], ...
    'LineWidth', 1.2);
without_line = plot([1 1] * marks.without, y_limits, ':', 'Color', [0.3 0.3 0.3], ...
    'LineWidth', 2);
text(marks.fitted + 0.01, label_y, sprintf('fitted %.2f', marks.fitted), ...
    'Color', [0.3 0.3 0.3], 'FontSize', 10);
text(marks.without + 0.01, label_y, sprintf('without %s %.2f', marks.without_mouse, ...
    marks.without), 'Color', [0.3 0.3 0.3], 'FontSize', 10);
handles = [fitted_line, band, without_line];
end

function draw_free_side(y_values, voxels, value_format, y_limits, is_log, colours)
% The region's two clusters on the index, beside a panel's axis and on its
% scale: the experimental-higher one filled, the control-higher one open, each
% with its value, or 'none' where the sign has no cluster.

hold on;
box on;
grid on;
point_colours = {colours.exp, colours.ctrl};
is_filled = [true, false];
for k = 1:2
    if is_filled(k)
        scatter(0.3, y_values(k), 60, point_colours{k}, 'filled');
    else
        scatter(0.3, y_values(k), 60, point_colours{k}, 'LineWidth', 1.4);
    end
    if voxels(k) == 0
        value_text = 'none';
    else
        value_text = sprintf(value_format, y_values(k));
    end
    text(0.45, y_values(k), value_text, 'Color', point_colours{k}, 'FontSize', 10);
end
if is_log
    set(gca, 'YScale', 'log');
    yticks([0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]);
    yline(0.05, '--', 'Color', [0.4 0.4 0.4]);
end
set(gca, 'FontSize', 11, 'YTickLabel', [], 'XTick', []);
ylim(y_limits);
xlim([0 1]);
xlabel({'no scale', '|L - R| / (L + R)'}, 'FontSize', 11);
end
