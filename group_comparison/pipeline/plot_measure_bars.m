function fig_bars = plot_measure_bars(perm, T_regions, bar_measure, exp_type, ...
    file_tag, map_titles, notes, a_priori_test)
%PLOT_MEASURE_BARS  The bars of one region measure, with the p of its permutation test.
%   fig_bars = PLOT_MEASURE_BARS(perm, T_regions, bar_measure, exp_type,
%   file_tag, map_titles, notes) draws the bars of one region measure of
%   region_permutation_test (perm), one panel per map, side by side, titled
%   by map_titles: the regions with a score other than zero, the largest
%   |score| at the top, signed (experimental higher to the right), grey by the
%   corrected permutation p, darker for a smaller p, starred below 0.05, a
%   dagger at an uncorrected p below 0.05, and the uncorrected p of the
%   regions named in advance (T_regions.a_priori) in bold beside their bars;
%   the regions expected to change after exp_type in bold magenta. The title
%   says which p tests what; notes, a cell of lines (empty for none), go
%   under it. Returns the figure, which the caller saves. Used by
%   group_differences and scale_free_test.
%
%   fig_bars = PLOT_MEASURE_BARS(..., a_priori_test) gives the regions named
%   in advance the p of a test other than the score's of either sign, such as
%   a one-sided one: a struct with p (regions x maps, NaN for a region not
%   named in advance) and label (what that p is, e.g. 'rws higher,
%   one-sided'). Their bold p and their dagger are then that p, with its label.

% the regions named in advance carry the uncorrected p of either sign, unless
% the test named in advance is another
if nargin < 8
    a_priori_test = [];
end

k = find(strcmp(perm.measure_names, bar_measure));
n_splits = size(perm.splits.in_ctrl, 1);
n_maps = numel(perm.maps);
panel_axes = gobjects(1, n_maps);
mark_texts = cell(1, n_maps);
n_drawn = zeros(1, n_maps);

fig_bars = figure('Visible', 'off', 'Name', ...
    ['Region_Surprise_Bar_' bar_measure '_' file_tag], 'Color', 'w', ...
    'Units', 'Normalized', 'Position', [0 0 0.9 0.9]);
for m = 1:n_maps
    score = perm.maps{m}.score(:, k);
    p_perm = perm.maps{m}.p_perm(:, k);
    p_fwer = perm.maps{m}.p_fwer(:, k);

    % the p of the regions named in advance: their uncorrected p, or that of the
    % test named in advance when it is another
    if isempty(a_priori_test)
        p_named = p_perm;
        named_label = '';
    else
        p_named = a_priori_test.p(:, m);
        named_label = a_priori_test.label;
    end

    % the regions with a score, the largest |score| last, which barh draws at the top
    drawn = find(~isnan(score) & score ~= 0);
    [~, order] = sort(abs(score(drawn)), 'ascend');
    drawn = drawn(order);

    % the bars, the label of the measure, and the regions expected to change in
    % bold magenta
    subplot(1, n_maps, m);
    mark_texts{m} = draw_measure_bars(score(drawn), p_perm(drawn), p_fwer(drawn), ...
        T_regions.label(drawn), T_regions.a_priori(drawn), n_splits, p_named(drawn), ...
        named_label);
    xlabel(measure_label(bar_measure, perm));
    title(sprintf('%s: %d of %d regions with a score', map_titles{m}, numel(drawn), ...
        height(T_regions)));
    highlight_surprise_regions(exp_type);
    panel_axes(m) = gca;
    n_drawn(m) = numel(drawn);
end

% the title, and what the marks and the magenta names mean: the corrected p tests
% a search over all the regions, the uncorrected p only a region named in advance
dagger = char(8224);
title_line = sprintf('Region %s, shaded by the corrected p of %d label permutations - %s', ...
    bar_measure, n_splits, strrep(file_tag, '_', ' '));
a_priori_acronyms = T_regions.acronym(T_regions.a_priori);
other_notes = {};
if isempty(a_priori_acronyms)
    marks_line = sprintf(['* corrected p < 0.05: the test of a search over all %d ' ...
                          'regions;   %s uncorrected p < 0.05: not a test, no region ' ...
                          'was named in advance'], height(T_regions), dagger);
elseif isempty(a_priori_test)
    marks_line = sprintf(['* corrected p < 0.05: the test of a search over all %d ' ...
                          'regions;   %s uncorrected p < 0.05: the test of a region ' ...
                          'named in advance, %s (its p in bold beside its bar)'], ...
                          height(T_regions), dagger, strjoin(a_priori_acronyms, ', '));
    other_notes{end + 1} = sprintf('a %s on any other region is not a test', dagger);
else
    marks_line = sprintf(['* corrected p < 0.05: the test of a search over all %d ' ...
                          'regions;   %s uncorrected p < 0.05: the test of a region ' ...
                          'named in advance, %s (its p, %s, in bold beside its bar)'], ...
                          height(T_regions), dagger, strjoin(a_priori_acronyms, ', '), ...
                          a_priori_test.label);
    other_notes{end + 1} = sprintf(['a %s on any other region is not a test (its p ' ...
                                    'of either sign)'], dagger);
end
if ~isempty(expected_regions(exp_type))
    other_notes{end + 1} = 'names in magenta: the regions expected to change';
end
title_lines = {title_line, ['\rm\fontsize{11}' marks_line]};
if ~isempty(other_notes)
    title_lines{end + 1} = ['\rm\fontsize{11}' strjoin(other_notes, ';   ')];
end
for i = 1:numel(notes)
    title_lines{end + 1} = ['\rm\fontsize{11}' notes{i}]; %#ok<AGROW>
end
sgtitle(title_lines, 'FontSize', 14, 'FontWeight', 'bold');

% in each panel, once the titles have taken their room: the daggers no taller
% than a row, the axis widened for the marks and the p beside the bars, and its
% ticks fixed
for m = 1:n_maps
    if n_drawn(m) > 0
        fit_daggers_to_rows(panel_axes(m), mark_texts{m}, n_drawn(m));
        widen_for_texts(panel_axes(m), mark_texts{m});
        set_value_ticks(panel_axes(m));
    end
end
end

% ===== Local functions =====

function mark_texts = draw_measure_bars(values, p_perm, p_fwer, labels, ...
    is_a_priori, n_splits, p_named, named_label)
% One panel of signed bars, grey by their corrected p, starred below 0.05, a
% dagger where only the uncorrected p is below 0.05, the uncorrected p in bold
% beside the bars of the regions named in advance, with the colour bar of the
% corrected p; returns the texts beside the bars, whose font size is the
% dagger's (the star and the p have their own). For the regions named in
% advance the uncorrected p is p_named, that of the test named in advance,
% written with named_label when it has one.

c_map = sep_palette('bars');
mark_texts = gobjects(0);
if isempty(values)
    text(0.5, 0.5, 'no region with a score', 'HorizontalAlignment', 'center');
    axis off;
    return
end

% the bars, each grey by its corrected p, the values written out on the axis (a
% common exponent would sit on the axis label)
b = barh(values, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = c_map(p_shade_index(p_fwer, n_splits, size(c_map, 1)), :);
ax = gca;
ax.XAxis.Exponent = 0;

% at the end of each bar, outside it: a star at a corrected p < 0.05, else a
% dagger at an uncorrected p < 0.05; for a region named in advance, its
% uncorrected p in bold beside the mark (the braces keep each size to its part)
hold on;
p_uncorrected = p_perm;
p_uncorrected(logical(is_a_priori)) = p_named(logical(is_a_priori));
for i = 1:numel(values)
    mark = '';
    if p_fwer(i) < 0.05
        mark = '{\fontsize{12}*}';
    elseif p_uncorrected(i) < 0.05
        mark = char(8224);
    end
    note = '';
    if is_a_priori(i) && isempty(named_label)
        note = sprintf('{\\fontsize{10}\\bf a priori: p = %.2g}', p_uncorrected(i));
    elseif is_a_priori(i)
        note = sprintf('{\\fontsize{10}\\bf a priori, %s: p = %.2g}', named_label, ...
            p_uncorrected(i));
    end
    if isempty(mark) && isempty(note)
        continue
    end

    % the mark next to the bar's end, the p after it
    if values(i) > 0
        beside = [' ' mark note];
        alignment = 'left';
    else
        beside = [note ' ' mark ' '];
        alignment = 'right';
    end
    mark_texts(end + 1) = text(values(i), i, beside, 'FontSize', 12, ...
        'HorizontalAlignment', alignment); %#ok<AGROW>
end

% the region names as tick labels, in the bars' order
yticks(1:numel(labels));
yticklabels(labels);
ylim([0 numel(labels) + 1]);
grid on;
set(gca, 'FontSize', 10);

% the colour bar of the corrected p
add_p_colorbar(c_map, n_splits);
end

function fit_daggers_to_rows(ax, texts, n_rows)
% Sets the texts beside the bars, whose size is the dagger's, to 12 points or
% less, so that daggers on adjacent rows stay apart.

% the dagger is about 0.9 of its font size tall, and the export at 300 dpi draws
% text about 10% larger than the screen: 0.8 of a row keeps a gap between rows
% (59 regions give rows of 9.6 points and 7.7-point daggers; at 12 points the
% daggers of adjacent rows join into one line)
drawnow;
ax.Units = 'points';
row_pt = ax.Position(4) / (n_rows + 1);
ax.Units = 'normalized';
set(texts, 'FontSize', min(12, 0.8 * row_pt));
end

function widen_for_texts(ax, texts)
% Widens the x axis of ax so that each text, anchored at the end of a bar, ends
% inside it, with 5% of the axis to spare.

% three passes, each on the figure as drawn: the new tick labels can take some
% of the axis's width
for pass = 1:3
    drawnow;
    for i = 1:numel(texts)
        limits = xlim(ax);

        % a text keeps its width on the screen, so it takes the same fraction of
        % the axis whatever the limits: the limit on its side moves until the rest
        % of the axis, room, holds the span from the other limit to the bar's end;
        % the width taken 15% larger, since the export at 300 dpi draws text up to
        % about 10% wider than the screen measures it
        width_fraction = 1.15 * texts(i).Extent(3) / diff(limits) + 0.05;
        room = 1 - width_fraction;
        anchor = texts(i).Position(1);
        if strcmp(texts(i).HorizontalAlignment, 'left')
            limits(2) = max(limits(2), limits(1) + (anchor - limits(1)) / room);
        else
            limits(1) = min(limits(1), limits(2) - (limits(2) - anchor) / room);
        end
        xlim(ax, limits);
    end
end
end

function set_value_ticks(ax)
% The x ticks of ax as MATLAB places them, but every other one, zero kept, when
% there are more than seven, which it would draw rotated; their labels as MATLAB
% writes them, but zero as 0, which it writes 0x10^0 when the others are powers
% of ten.

tick_values = xticks(ax);
if numel(tick_values) > 7
    step = tick_values(2) - tick_values(1);
    is_even = mod(round(tick_values / step), 2) == 0;
    tick_values = tick_values(is_even);
end
xticks(ax, tick_values);
tick_labels = xticklabels(ax);
tick_labels(tick_values == 0) = {'0'};
xticklabels(ax, tick_labels);
end

function label = measure_label(measure, perm)
% The axis label of a region measure.

settings = perm.settings;
switch measure
    case 'share'
        label = sprintf('fraction of the voxels with a t at p < %g', settings.p_thresh);
    case 'sum'
        label = sprintf('summed surprise (-log_{10} p) of the voxels at p < %g, signed', ...
            settings.p_thresh);
    case 'topvol'
        label = sprintf('mean surprise of the top %g mm^3 (%d voxels), signed', ...
            settings.topvol_mm3, perm.topvol_k);
    case 'cluster'
        label = sprintf(['mass (summed surprise) of the heaviest cluster at p < %g, ' ...
                         'signed'], settings.cluster_p);
    otherwise
        label = sprintf('%g quantile of the surprise, signed', settings.region_quantile);
end
end
