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
%                                    either sign
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
%   step 2. Two panels against it, at the slopes of the sweep and the fitted
%   one: the one-sided p of the heaviest cluster of each sign over every split
%   of the mice (for the experimental-higher sign the test the scale-free test
%   named in advance), and each cluster's voxels. A slope of the sweep within
%   half its step of the fitted one is left out, so no two markers sit on each
%   other. Both panels mark the fitted slope with its jackknife standard error,
%   from the slopes s_i refitted without each of the n mice,
%       se = sqrt((n - 1) / n * sum((s_i - mean(s)).^2))
%   and the slope refitted without the mouse that moves it most. Under the
%   title, the experimental-higher cluster's p at the fitted slope, one-sided
%   and of either sign (step 3's p, the one quoted), and its one-sided p on
%   the index |L - R| / (L + R) of run_scale_free_test, where no scale enters.
%   Every number is read from the tables. The fitted slope must be fold 0's,
%   and its cluster step 3's of the scale-free test's table, so the two steps
%   read the same step 3.

% settings of run_rescaling_figure, under the names the code below uses
exp_type = run_settings.exp_type;
ctrl_name = run_settings.ctrl_name;
exp_name = run_settings.exp_name;
region_name = run_settings.region_name;
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

%% The slopes drawn and their marks

% the sweep's slopes and the fitted one, less the sweep's slope beside the
% fitted one (1.5 beside 1.49)
drawn = drawn_slopes(T_slopes);

% the slope refitted without each mouse, its jackknife standard error, and the
% mouse whose leaving out moves it most
marks = slope_marks(T_folds, fitted.slope);

% the region's experimental-higher cluster on the index, where no scale enters
free_exp = cluster_row(T_free, [exp_type ' higher, index'], free_file);

%% Figure

names = struct('ctrl', ctrl_name, 'exp', exp_name, 'region', region_name);
result = result_line(fitted, free_exp, names);
print_summary(drawn, marks, result, names);
plot_figure(drawn, marks, result, names, file_tag, comp_out_dir);
fprintf('Rescaling figure saved to: %s\n', ...
    fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.png']));

end

% ===== Local functions: the tables and their numbers =====

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

function drawn = drawn_slopes(T_slopes)
% The rows drawn, in order of the slope: the fitted one and the sweep's, less
% any of the sweep's within half the sweep's step of the fitted one.

is_sweep = T_slopes.is_fitted == 0;
sweep_step = min(diff(T_slopes.slope(is_sweep)));
fitted_slope = T_slopes.slope(T_slopes.is_fitted == 1);
is_beside = is_sweep & abs(T_slopes.slope - fitted_slope) < sweep_step / 2;
drawn = T_slopes(~is_beside, :);
end

function marks = slope_marks(T_folds, fitted_slope)
% The fitted slope; the slopes refitted without each mouse and their jackknife
% standard error; the mouse whose leaving out moves the slope most, and its
% slope.

is_fold = T_folds.fold > 0;
refitted = T_folds.slope(is_fold);
left_out = T_folds.left_out(is_fold);
n_mice = numel(refitted);
[~, i_most] = max(abs(refitted - fitted_slope));

marks = struct();
marks.fitted = fitted_slope;
marks.se = sqrt((n_mice - 1) / n_mice * sum((refitted - mean(refitted)).^2));
marks.n_mice = n_mice;
marks.without = refitted(i_most);
marks.without_mouse = strrep(left_out{i_most}, '_Gria1', '');
end

function result = result_line(fitted, free_exp, names)
% The line under the title: the experimental-higher cluster's p at the fitted
% slope, one-sided and of either sign, and on the index.

result = sprintf(['At the fitted slope %.2f the %s-higher patch has p %.3f one-sided ' ...
                  '(%.3f for either sign, the value quoted). On the ratio ' ...
                  '|L - R| / (L + R), where no rescaling enters, it has p %.2f.'], ...
                  fitted.slope, names.exp, fitted.p_positive, fitted.p_either_sign, ...
                  free_exp.p_one_sided);
end

function print_summary(drawn, marks, result, names)
% The numbers the figure shows.

fprintf('%s, %s against %s, the heaviest cluster of each sign over %d splits:\n', ...
    names.region, names.exp, names.ctrl, drawn.n_splits(1));
for s = 1:height(drawn)
    note = '';
    if drawn.is_fitted(s) == 1
        note = ' (fitted)';
    end
    fprintf(['  slope %.2f%s: %s higher %d voxels, p %.3f; %s higher %d voxels, ' ...
             'p %.3f\n'], drawn.slope(s), note, names.exp, drawn.cluster_n(s), ...
             drawn.p_positive(s), names.ctrl, drawn.negative_n(s), drawn.p_negative(s));
end
fprintf(['  fitted slope %.4f, jackknife standard error %.3f over %d mice; refitted ' ...
         'without %s %.4f\n'], marks.fitted, marks.se, marks.n_mice, ...
         marks.without_mouse, marks.without);
fprintf('  %s\n', result);
end

% ===== Local functions: figure =====

function plot_figure(drawn, marks, result, names, file_tag, comp_out_dir)
% The p and the size of each sign's cluster against the slope, the title with
% the result under it, and one legend under both panels.

% each sign's colour: red where the experimental group is the higher, dark
% grey where the control group is
colours = struct();
colours.exp = sep_palette('experimental_higher');
colours.ctrl = sep_palette('control_higher');
x_limits = [min(drawn.slope) - 0.05, max(drawn.slope) + 0.05];

% a fixed size in inches, so the fonts keep their size against the panels on
% any screen
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0 0 16 7]);

axes('Position', [0.065 0.2 0.475 0.63]);
handles = draw_p_panel(drawn, marks, colours, x_limits, names);
axes('Position', [0.625 0.2 0.36 0.63]);
draw_size_panel(drawn, marks, colours, x_limits, names);

% one legend, under both panels; no-break spaces after the first name keep it
% apart from the second's line (the legend drops plain spaces at a name's end)
labels = {
    sprintf('%s-higher patch (one-sided p, the test named in advance)%s', names.exp, ...
        repmat(char(160), 1, 8))
    sprintf('%s-higher patch', names.ctrl)
    };
lgd = legend(handles, labels, 'Orientation', 'horizontal', 'FontSize', 13, ...
    'Box', 'off');
lgd.Units = 'normalized';
lgd.Position(1) = 0.5 - lgd.Position(3) / 2;
lgd.Position(2) = 0.015;

title_line = sprintf(['%s against %s, %s: the p of the test named in advance ' ...
                      'depends on the rescaling'], names.exp, names.ctrl, names.region);
sgtitle({title_line, ['\fontsize{12}' result]}, 'FontSize', 16, ...
    'FontWeight', 'normal');
saveas(fig, fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Rescaling_Test_' file_tag '.png']), ...
    'Resolution', 300);
end

function handles = draw_p_panel(drawn, marks, colours, x_limits, names)
% Each sign's p against the slope on a log axis, the slope's marks behind it
% and 0.05 dotted; the marks named at the panel's foot. Returns the two series
% the legend names.

p_limits = [2e-3 2];
hold on;
draw_slope_marks(marks, x_limits, p_limits);
yline(0.05, ':', '0.05', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
    'FontSize', 12, 'LabelHorizontalAlignment', 'left');
exp_line = draw_series(drawn.slope, drawn.p_positive, 'o', colours.exp);
ctrl_line = draw_series(drawn.slope, drawn.p_negative, 's', colours.ctrl);

% each mark named right of its line, at the foot, under the curves
label_y = 1.2 * p_limits(1);
text(marks.fitted + 0.012, label_y, {sprintf('fitted slope %.2f', marks.fitted), ...
    sprintf('shading: \\pm1 SE (%.2f)', marks.se)}, 'VerticalAlignment', 'bottom', ...
    'FontSize', 12);
text(marks.without + 0.012, label_y, {'without', marks.without_mouse}, ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', [0.4 0.4 0.4]);

set(gca, 'YScale', 'log');
ylim(p_limits);
yticks([0.01 0.1 1]);
yticklabels({'0.01', '0.1', '1'});
finish_panel(x_limits, names, 'A. p against the rescaling');
ylabel(sprintf('cluster-mass p (%d relabellings)', drawn.n_splits(1)));
handles = [exp_line, ctrl_line];
end

function draw_size_panel(drawn, marks, colours, x_limits, names)
% Each sign's cluster voxels against the slope, the slope's marks behind them.

size_limits = [0, 1.05 * max([drawn.cluster_n; drawn.negative_n])];
hold on;
draw_slope_marks(marks, x_limits, size_limits);
draw_series(drawn.slope, drawn.cluster_n, 'o', colours.exp);
draw_series(drawn.slope, drawn.negative_n, 's', colours.ctrl);

ylim(size_limits);
finish_panel(x_limits, names, 'B. which patch is larger');
ylabel(sprintf('patch size in the %s (voxels)', names.region));
end

function draw_slope_marks(marks, x_limits, y_limits)
% The fitted slope in black over a light grey band of +/- 1 jackknife standard
% error, cut at the axis, and the slope refitted without the mouse that moves
% it most, grey and dashed.

band_x = [max(marks.fitted - marks.se, x_limits(1)), ...
    min(marks.fitted + marks.se, x_limits(2))];
patch(band_x([1 2 2 1]), y_limits([1 1 2 2]), [0.93 0.93 0.93], ...
    'EdgeColor', 'none');
plot([1 1] * marks.fitted, y_limits, '-', 'Color', [0 0 0], 'LineWidth', 1.2);
plot([1 1] * marks.without, y_limits, '--', 'Color', [0.6 0.6 0.6], ...
    'LineWidth', 1.5);
end

function handle = draw_series(x, y, marker, colour)
% One series, its markers filled and joined by a line.

handle = plot(x, y, ['-' marker], 'Color', colour, 'MarkerFaceColor', colour, ...
    'LineWidth', 2, 'MarkerSize', 8);
end

function finish_panel(x_limits, names, panel_title)
% The slope axis, the panel's title on the left, open axes with ticks out.

set(gca, 'FontSize', 13, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top', ...
    'TitleHorizontalAlignment', 'left');
xlim(x_limits);
xtickformat('%.1f');
xlabel(sprintf('slope of the line aligning %s onto %s', names.exp, names.ctrl));
title(panel_title, 'FontWeight', 'normal');
end
