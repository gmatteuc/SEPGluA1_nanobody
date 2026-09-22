function key = cohort_atlas_key(group, ages)
%COHORT_ATLAS_KEY Which atlas a cohort is registered to.
%
%   key = cohort_atlas_key('naive', [])      -> 'ccf'
%   key = cohort_atlas_key('young', 20)      -> 'demba_p20'
%
% The decision itself is recorded in get_atlas and in the project notes:
% adults stay on the Allen CCF, the young brains go to the DeMBA atlas of
% their own age. An age whose atlas folder has not been built is an error
% here rather than a silent fall-back onto a neighbouring age or an adult
% atlas -- registering a P16 brain to a P20 (or adult) template would be a
% quiet way to manufacture a developmental difference.

if nargin < 2, ages = []; end

switch lower(group)
    case {'naive', 'rws', 'behavior'}
        key = 'ccf';
    case 'young'
        if isempty(ages)
            error(['cohort_atlas_key: the young group spans several ages, each on ' ...
                   'its own atlas. Give an age, e.g. ''young_P20''.']);
        end
        u = unique(ages(:))';
        if numel(u) > 1
            error(['cohort_atlas_key: ages P%s would need one atlas each. ' ...
                   'Analyse one age at a time.'], ...
                   strjoin(arrayfun(@num2str, u, 'UniformOutput', false), ' and P'));
        end
        key = sprintf('demba_p%d', u);
        p = get_paths();
        if ~exist(fullfile(p.data, sprintf('atlas_demba_p%d', u)), 'dir')
            error(['cohort_atlas_key: no DeMBA atlas has been built for P%d yet.\n' ...
                   'Build it:  tools\\venv_atlas\\Scripts\\python.exe build_demba_atlas.py %d'], u, u);
        end
    otherwise
        error('cohort_atlas_key: unknown group ''%s''.', group);
end
end
