function p = get_paths()
%GET_PATHS  Where everything in the project lives.
%   P = GET_PATHS() returns the project's folders, worked out from where this
%   file sits rather than written out, so the whole tree can be moved to
%   another drive, or handed over on an external disk, and still run without
%   anyone editing a path. The only literals are the production folders of
%   the guard below.
%
%   Fields of P:
%     root             project root, the parent of this code folder
%     data             derived data, one folder per group
%     atlas            Allen atlas volumes and ontology CSVs
%     lightsuite       vendored LightSuite (third_party\LightSuite)
%     yaml             vendored yamlmatlab (third_party\yamlmatlab)
%     elastix          vendored matlab_elastix (third_party\matlab_elastix)
%     bioformats       vendored Bio-Formats reader (third_party\BioformatsImage)
%     auto_annotation  the automatic annotation's Python engine
%                      (registration\auto_annotation)
%     code             this folder
%
%   The layout is a plain sibling arrangement:
%     <root>\code\    this file, and the code in its pipeline folders
%     <root>\data\    everything the pipeline produces
%
%   The environment variable SEP_DATA_ROOT points the whole data tree
%   elsewhere. The refactor's checks run on copies of the data and must never
%   write into the real one (docs/REFACTOR_PLAN.md), so a copy of the code is
%   refused the production data, and neither the snapshot on G: nor the code
%   folder is accepted as a data root. The Python route applies the same rules
%   (mapping/sepmap/config.py, and atlas/build_demba_atlas.py, which sits
%   outside the package).
%
%   Raw .czi files are copied into <root>\data\<group>\<mouse>\ before
%   processing, so the read-only lab share is not part of this at all (see
%   run_copy_raw_data).

p.code = fileparts(mfilename('fullpath'));
p.root = fileparts(p.code);

% the data root, moved by SEP_DATA_ROOT and refused where a check must not write
p.data = fullfile(p.root, 'data');
if ~isempty(getenv('SEP_DATA_ROOT'))
    p.data = canonical(getenv('SEP_DATA_ROOT'));
end
check_data_root(p.code, p.data);
p.atlas = fullfile(p.data, 'atlas');

% the vendored toolboxes and the automatic annotation's engine, in the code folder
p.lightsuite = fullfile(p.code, 'third_party', 'LightSuite');
p.yaml       = fullfile(p.code, 'third_party', 'yamlmatlab');
p.elastix    = fullfile(p.code, 'third_party', 'matlab_elastix');
p.bioformats = fullfile(p.code, 'third_party', 'BioformatsImage');
p.auto_annotation = fullfile(p.code, 'registration', 'auto_annotation');

end

% ===== Local functions =====

function check_data_root(code_dir, data_dir)
% Refuse a data root inside the production data (for a copy of the code), the
% G: snapshot or the code folder.

% a copy of the code placed next to the production code folder would work out
% the production data as its own, and a check run from it would overwrite real outputs
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

% the snapshot on G: is a backup, never a data root
if startsWith(lower(data_dir), lower('G:\sep_histology_snapshot'))
    error('get_paths: the data root %s is inside the snapshot on G:, which is a backup.', ...
        data_dir);
end

% a copy of the code in a worktree inside the code folder would otherwise put
% its data, and the folders drivers create, inside the code tree
if strcmpi(data_dir, production_code) || ...
        startsWith(lower(data_dir), lower([production_code '\']))
    error(['get_paths: the data root %s is inside the code folder (a copy of the code in ' ...
           'a worktree there?). Set SEP_DATA_ROOT to its check tree.'], data_dir);
end

end

function s = canonical(s)
% One spelling for every path, so the comparisons above cannot be dodged: absolute,
% '.', '..' and doubled separators resolved, backslashes, no trailing separator.

s = strrep(char(s), '/', '\');
if isempty(regexp(s, '^([A-Za-z]:\\|\\\\)', 'once'))
    s = fullfile(pwd, s);
end
s = char(java.io.File(s).toPath().normalize().toString());
while numel(s) > 3 && s(end) == '\'
    s(end) = [];
end

end
