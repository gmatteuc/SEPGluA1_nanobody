function shade = p_shade_index(p, n_splits, n_levels)
%P_SHADE_INDEX  The grey level of a corrected p, for the region measures' figures.
%   shade = P_SHADE_INDEX(p, n_splits, n_levels) returns the grey level of each
%   corrected p: 1 (the lightest) at p = 1 to n_levels (black) at the smallest
%   p the splits allow, 1 / n_splits, on a log scale. add_p_colorbar draws its
%   colour bar. Used by plot_measure_bars and group_differences.

darkness = log(p) / log(1 / n_splits);
darkness(isnan(darkness)) = 0;
darkness = min(max(darkness, 0), 1);
shade = max(1, ceil(darkness * n_levels));
end
