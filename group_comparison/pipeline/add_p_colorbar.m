function add_p_colorbar(c_map, n_splits)
%ADD_P_COLORBAR  The colour bar of the corrected p, for the region measures' figures.
%   ADD_P_COLORBAR(c_map, n_splits) gives the current axes the colormap c_map
%   and its colour bar of the corrected p, from 1 (light, at the bottom) to
%   1 / n_splits (black, at the top), its ticks where p_shade_index shades
%   them. Used by plot_measure_bars and group_differences.

colormap(gca, c_map);
clim([0 1]);
cb = colorbar;

% the ticks at these p, placed as p_shade_index shades them, from the bottom
tick_p = [1, 0.5, 0.1, 0.05, 0.01, 1 / n_splits];
tick_p = sort(unique(tick_p(tick_p >= 1 / n_splits)), 'descend');
cb.Ticks = log(tick_p) / log(1 / n_splits);
cb.TickLabels = arrayfun(@(p) sprintf('%.3g', p), tick_p, 'UniformOutput', false);
cb.Label.String = 'corrected p';
end
