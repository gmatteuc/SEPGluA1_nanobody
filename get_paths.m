function p = get_paths()
% GET_PATHS Where everything in the project lives.
%
%   p.root        project root, the parent of this code folder
%   p.data        derived data, one folder per group
%   p.atlas       Allen atlas volumes and ontology CSVs
%   p.lightsuite  vendored LightSuite
%   p.yaml        vendored yamlmatlab
%   p.elastix     vendored matlab_elastix
%   p.bioformats  vendored Bio-Formats reader
%   p.code        this folder
%
% The layout is a plain sibling arrangement:
%
%   <root>\code\    this file, and the pipeline scripts
%   <root>\data\    everything the pipeline produces
%
% All of it is worked out from where this file sits rather than written out as
% a literal. That way the whole tree can be moved to another drive, or handed
% to someone else on an external disk, and still run without anyone editing
% paths. Nothing here says D: or E:.
%
% Raw .czi are copied into <root>\data\<group>\<mouse>\ before processing, so
% the read-only lab share is not part of this at all -- see P0_copy_raw_data.

p.code = fileparts(mfilename('fullpath'));
p.root = fileparts(p.code);

p.data  = fullfile(p.root, 'data');
p.atlas = fullfile(p.data, 'atlas');

p.lightsuite = fullfile(p.code, 'LightSuite-main');
p.yaml       = fullfile(p.code, 'yamlmatlab');
p.elastix    = fullfile(p.code, 'matlab_elastix-master');
p.bioformats = fullfile(p.code, 'BioformatsImage');

end
