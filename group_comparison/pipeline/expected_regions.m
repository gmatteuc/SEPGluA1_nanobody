function names = expected_regions(exp_type)
%EXPECTED_REGIONS  The regions expected to change in the plasticity comparison.
%   names = EXPECTED_REGIONS(exp_type) returns the names of the regions
%   expected to change after exp_type, which the bars label in bold magenta:
%   the whisker system, as in the coarse bars, and the rostrolateral visual
%   area; none for another group. Used by highlight_surprise_regions,
%   plot_measure_bars and group_differences.

switch exp_type
    case {'rws', 'behavior'}
        names = {'Primary somatosensory area, barrel field', ...
            'Ventral posteromedial nucleus of the thalamus', ...
            'Posterior complex of the thalamus', ...
            'Supplemental somatosensory area', 'Zona incerta', ...
            'Rostrolateral visual area'};
    otherwise
        names = {};
end
end
