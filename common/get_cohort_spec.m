function S = get_cohort_spec(spec)
%GET_COHORT_SPEC Resolve a cohort name into everything the analysis needs.
%
%   S = get_cohort_spec('naive')        an adult group, as always
%   S = get_cohort_spec('young_P20')    the young group, one age only
%   S = get_cohort_spec('young_P16P20') the young group, several ages
%
% Returns a struct:
%   group      registry group name ('naive', 'rws', 'behavior', 'young')
%   ages       age filter in days ([] for a whole group)
%   tag        filename suffix P5 uses for an age-filtered aggregate: '' for
%              a whole group, '_P20' for young_P20 -- so nano_4d<tag>.mat,
%              collected_mice<tag>.mat, nano_4d_normalized<tag>.mat
%   base_dir   <data>\<group>
%   atlas_key  'ccf' for the adults, 'demba_p20' for P20 (see cohort_atlas_key)
%   mice       the mouse names, in the order P5 stacked them: from
%              collected_mice<tag>.mat when P5 wrote one, else from the
%              registry (the adults were aggregated before P5 recorded that)
%   label      the spec string, for titles and folder names
%
% This is the one place that knows how a cohort spec maps onto files, so the
% analysis scripts (P6bis, P8) can take 'young_P20' and 'naive' alike and
% never open an atlas or list mice themselves. An adult spec resolves to
% exactly the names and files those scripts used before.

spec = char(spec);
tok = regexp(spec, '^(\w+?)(_P[\dP]+)?$', 'tokens', 'once');
if isempty(tok)
    error('get_cohort_spec: cannot parse cohort spec ''%s''.', spec);
end
group = tok{1};
tag   = tok{2};
if isempty(tag)
    ages = [];
else
    ages = cellfun(@str2double, regexp(tag, '\d+', 'match'));
end

paths  = get_paths();
cohort = get_cohort('groups', {group});
if ~isempty(ages)
    cohort = cohort(ismember([cohort.age_days], ages));
end
if isempty(cohort)
    error('get_cohort_spec: no mice match ''%s''.', spec);
end

S = struct();
S.group     = group;
S.ages      = ages;
S.tag       = tag;
S.label     = spec;
S.base_dir  = fullfile(paths.data, group);
S.atlas_key = cohort_atlas_key(group, ages);

collected = fullfile(S.base_dir, ['collected_mice' tag '.mat']);
if exist(collected, 'file')
    C = load(collected, 'collected_mice');
    S.mice = C.collected_mice;
else
    S.mice = {cohort.name};
end
end
