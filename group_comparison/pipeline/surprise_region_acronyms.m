function acronyms = surprise_region_acronyms(T_members)
%SURPRISE_REGION_ACRONYMS  The regions of the surprise bars, by atlas acronym.
%   acronyms = SURPRISE_REGION_ACRONYMS(T_members) returns the isocortical
%   areas, taken from the atlas (T_members, the table of
%   parcellation_to_parcellation_term_membership.csv), then the declared
%   subcortical regions, as one column of acronyms. Used by surprise_regions
%   and by group_differences, which checks the regions named in advance.

acronyms = [isocortical_areas(T_members); surprise_subcortical_list()];
end

% ===== Local functions =====

function iso_acronyms = isocortical_areas(T_members)
% The acronyms of the atlas's structure-level terms in the isocortex, in the
% atlas's order: the cortical areas, each with its layers below it.

% the atlas values of the isocortex division
is_iso = strcmp(T_members.parcellation_term_set_name, 'division') & ...
    strcmp(T_members.parcellation_term_acronym, 'Isocortex');
iso_values = T_members.parcellation_index(is_iso);

% their terms at the structure level: one level for every area, so S1 comes as
% its seven subfields (SSp-n to SSp-un) and the anterior cingulate as ACAd and
% ACAv, where a deeper or shallower level of the ontology would mix the two
is_iso_area = strcmp(T_members.parcellation_term_set_name, 'structure') & ...
    ismember(T_members.parcellation_index, iso_values);
iso_acronyms = unique(T_members.parcellation_term_acronym(is_iso_area), 'stable');
end

function sub_acronyms = surprise_subcortical_list()
% The regions of the surprise bars outside the isocortex, by atlas acronym,
% grouped by their atlas division in the atlas's order.

% olfactory areas: piriform cortex
olfactory = {'PIR'};

% hippocampal formation, without the subiculum, which has its own bar
hippocampal = {'HPF'; 'SUB'};

% cortical subplate: claustrum, basolateral amygdala
subplate = {'CLA'; 'BLA'};

% striatum: olfactory tubercle, nucleus accumbens, caudoputamen
striatum = {'OT'; 'ACB'; 'CP'};

% pallidum: external globus pallidus
pallidum = {'GPe'};

% thalamus: the whisker relays VPM and PO, the other sensory-motor and
% higher-order nuclei, the reticular nucleus and the dorsal geniculate group
thalamus = {'VPM'; 'VPL'; 'VM'; 'PO'; 'LP'; 'LD'; 'VAL'; 'MD'; 'PF'; 'RE'; 'CL'; ...
    'RT'; 'GENd'};

% hypothalamus, without the subthalamic nucleus and the zona incerta, which have
% their own bars
hypothalamus = {'HY'; 'STN'; 'ZI'};

% midbrain: its motor part, without the motor superior colliculus, which has its
% own bar, and the sensory superior colliculus
midbrain = {'MBmot'; 'SCm'; 'SCs'};

sub_acronyms = [olfactory; hippocampal; subplate; striatum; pallidum; thalamus; ...
    hypothalamus; midbrain];
end
