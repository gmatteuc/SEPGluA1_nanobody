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
% paths. The only literals are the production folders in the guard below.
%
% The environment variable SEP_DATA_ROOT points the whole data tree somewhere
% else. The refactor's checks run on copies of the data and must never write
% into the real one (docs/REFACTOR_PLAN.md), so a copy of the code is refused
% the production data, and the snapshot on G: is refused as a data root. The
% Python route applies the same rules (v2_paths.py).
%
% Raw .czi are copied into <root>\data\<group>\<mouse>\ before processing, so
% the read-only lab share is not part of this at all -- see P0_copy_raw_data.

p.code = fileparts(mfilename('fullpath'));
p.root = fileparts(p.code);

p.data = fullfile(p.root, 'data');
if ~isempty(getenv('SEP_DATA_ROOT'))
    p.data = canonical(getenv('SEP_DATA_ROOT'));
end
check_data_root(p.code, p.data);
p.atlas = fullfile(p.data, 'atlas');

p.lightsuite = fullfile(p.code, 'LightSuite-main');
p.yaml       = fullfile(p.code, 'yamlmatlab');
p.elastix    = fullfile(p.code, 'matlab_elastix-master');
p.bioformats = fullfile(p.code, 'BioformatsImage');

end


function check_data_root(code_dir, data_dir)
% A copy of the code placed next to the production code folder would work
% out the production data folder as its own, and a check run from it would
% overwrite real outputs; the snapshot on G: is a backup, never a data root.
production_code = 'D:\sep_histology\code';
production_data = 'D:\sep_histology\data';
code_dir = canonical(code_dir);
data_dir = canonical(data_dir);
inside_production = strcmpi(data_dir, production_data) || ...
    startsWith(lower(data_dir), lower([production_data '\']));
if inside_production && ~strcmpi(code_dir, production_code)
    error(['get_paths: this copy of the code (%s) would use the production data (%s). ' ...
           'Set SEP_DATA_ROOT to the data of its own check tree.'], code_dir, data_dir);
end
if startsWith(lower(data_dir), lower('G:\sep_histology_snapshot'))
    error('get_paths: the data root %s is inside the snapshot on G:, which is a backup.', data_dir);
end
% A copy of the code in a worktree inside the code folder would otherwise
% put its data, and the folders drivers create, inside the code tree.
if strcmpi(data_dir, production_code) || startsWith(lower(data_dir), lower([production_code '\']))
    error(['get_paths: the data root %s is inside the code folder (a copy of the code in ' ...
           'a worktree there?). Set SEP_DATA_ROOT to its check tree.'], data_dir);
end
end


function s = canonical(s)
% One spelling for every path, so the comparisons above cannot be dodged: a
% relative path is taken from the current folder, '.', '..' and doubled
% separators are resolved, backslashes, no trailing separator.
s = strrep(char(s), '/', '\');
if isempty(regexp(s, '^([A-Za-z]:\\|\\\\)', 'once'))
    s = fullfile(pwd, s);
end
s = char(java.io.File(s).toPath().normalize().toString());
while numel(s) > 3 && s(end) == '\'
    s(end) = [];
end
end
