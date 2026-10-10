function highlight_surprise_regions(exp_type)
%HIGHLIGHT_SURPRISE_REGIONS  The regions expected to change, in bold magenta on the bars.
%   HIGHLIGHT_SURPRISE_REGIONS(exp_type) recolours the tick labels of the
%   current axes: those of the regions expected to change after exp_type
%   (expected_regions) in bold magenta, the others in black. Used by
%   plot_measure_bars and group_differences.

% the regions expected to change
highlighted_areas = expected_regions(exp_type);

% each tick label recoloured with TeX markup, underscores as spaces
ax = gca;
ytl = ax.YTickLabel;
colored_labels = repmat({''}, size(ytl));
for i = 1:length(ytl)
    if ismember(ytl{i}, highlighted_areas)
        colored_labels{i} = ['\color{magenta} \bf ' strrep(ytl{i}, '_', ' ')];
    else
        colored_labels{i} = ['\color{black} ' strrep(ytl{i}, '_', ' ')];
    end
end
ax.YTickLabel = colored_labels;
end
