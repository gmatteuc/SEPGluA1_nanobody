function key = cohort_atlas_key(group, ages)
%COHORT_ATLAS_KEY Which atlas a cohort is registered to.
%
%   key = cohort_atlas_key('naive', [])      -> 'ccf'
%   key = cohort_atlas_key('young', 20)      -> 'demba_p20'
%
% The decision itself is recorded in get_atlas and in the project notes:
% adults stay on the Allen CCF, the young brains go to the DeMBA atlas of
% their own age. Only the P20 DeMBA has been built so far, so any other
% young age is an error here rather than a silent fall-back onto an adult
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
        if isequal(sort(ages(:))', 20)
            key = 'demba_p20';
        else
            error(['cohort_atlas_key: no DeMBA atlas has been built for P%s yet. ' ...
                   'Only demba_p20 exists (see get_atlas).'], ...
                   strjoin(arrayfun(@num2str, ages, 'UniformOutput', false), '/P'));
        end
    otherwise
        error('cohort_atlas_key: unknown group ''%s''.', group);
end
end
