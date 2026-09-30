function atlas = get_atlas(atlas_key)
% GET_ATLAS Resolve reference atlas location and parameters by key.
%
%   atlas = get_atlas()      returns the default ('ccf')
%   atlas = get_atlas('ccf') Allen Mouse Brain CCFv3, 10 um
%
%   Returned struct fields:
%     key              short identifier used in output tags
%     dir              folder holding the atlas volumes
%     template_file    grayscale average template (used for registration)
%     annotation_file  region-label volume (used for ROI masks)
%     boundary_file    region-boundary volume (used for overlays)
%     res_um           isotropic voxel size in micrometres
%     age_days         postnatal age this atlas represents (56 = adult). An
%                      age-matched atlas is only valid for mice of that age --
%                      a P36 brain does not belong on the P20 template.
%     default_aplims   AP crop used historically for this atlas
%     description      human-readable note
%
%   'ccf' reproduces exactly the paths and files every existing result was
%   produced with, so switching call sites from the hardcoded
%   'D:\sep_histology\data\atlas' literal to get_atlas('ccf') is a pure
%   refactor with no behavioural change.
%
%   NOTE ON default_aplims: the AP crop is currently read from each mouse's
%   local_settings.txt (`atlasaplims`), NOT from here. The value below is
%   recorded for reference and for future atlases only -- P1-P4 must keep
%   using sliceinfo.atlasaplims so existing behaviour is preserved.
%
%   ADDING A YOUNG-BRAIN ATLAS: register a new key here (e.g. 'devccf_p14')
%   pointing at its template/annotation volumes. Cross-group comparison then
%   requires both atlases to resolve into a common space -- DevCCF ships
%   CCF-linked labels for exactly this purpose. Do not swap atlases for a
%   comparison without that mapping.

if nargin < 1 || isempty(atlas_key)
    atlas_key = 'ccf';
end

p = get_paths();

switch lower(atlas_key)

    case 'ccf'
        atlas.key             = 'ccf';
        atlas.dir             = p.atlas;
        atlas.template_file   = 'average_template_10.nii.gz';
        atlas.annotation_file = 'annotation_10.nii.gz';
        atlas.boundary_file   = 'annotation_boundary_10.nii.gz';
        atlas.res_um          = 10;
        atlas.age_days        = 56;
        atlas.default_aplims  = [180 1079];
        atlas.description     = 'Allen Mouse Brain Common Coordinate Framework v3, 10 um (adult, P56)';

    otherwise
        % Any DeMBA age: 'demba_p16', 'demba_p20', 'demba_p28' ... One entry
        % rather than one case per age, because the folders are built to a
        % fixed recipe by build_demba_atlas.py and differ only in the age and
        % in the AP crop that age measures. An age whose folder has not been
        % built is an error, not a fall-back onto a neighbouring age.
        tok = regexp(lower(atlas_key), '^demba_p(\d+)$', 'tokens', 'once');
        if isempty(tok)
            error(['get_atlas: unknown atlas key "%s". Known keys: ''ccf'' and ' ...
                   '''demba_pN'' for any age N built by build_demba_atlas.py.'], atlas_key);
        end
        age = str2double(tok{1});
        atlas.key             = sprintf('demba_p%d', age);
        atlas.dir             = fullfile(p.data, sprintf('atlas_demba_p%d', age));
        atlas.template_file   = 'average_template_10.nii.gz';
        atlas.annotation_file = 'annotation_10.nii.gz';
        atlas.boundary_file   = '';
        atlas.res_um          = 20;
        atlas.age_days        = age;
        atlas.description     = sprintf(['DeMBA P%d (Carey 2025), Allen CCFv3 labels, ' ...
                                         '20 um isotropic'], age);

        if ~exist(atlas.dir, 'dir')
            error(['get_atlas: no atlas built for P%d.\n  %s does not exist.\n' ...
                   'Build it first:  tools\\venv_atlas\\Scripts\\python.exe build_demba_atlas.py %d'], ...
                   age, atlas.dir, age);
        end
        % The AP crop is measured per age when the folder is built (the two
        % methods are described in build_demba_atlas.py) and stored beside the
        % volumes, so it cannot drift away from the atlas it belongs to.
        aplims_file = fullfile(atlas.dir, 'aplims.txt');
        if ~exist(aplims_file, 'file')
            error(['get_atlas: %s is missing. Re-run build_demba_atlas.py %d, or write the ' ...
                   'two AP crop planes into that file.'], aplims_file, age);
        end
        lims = sscanf(fileread(aplims_file), '%d')';
        if numel(lims) ~= 2 || lims(2) <= lims(1)
            error('get_atlas: %s should hold two increasing plane numbers, found "%s".', ...
                aplims_file, strtrim(fileread(aplims_file)));
        end
        atlas.default_aplims  = lims;

        % Two things about these entries are deliberate and easy to get wrong:
        %
        % The files are named *_10.nii.gz but hold 20 um data. LightSuite finds
        % the atlas with which('average_template_10.nii.gz') in fourteen
        % different files, so the name is forced by the vendored code. The real
        % resolution is res_um here and px_atlas in local_settings.txt, and both
        % must say 20 for a young brain or the AP scale silently doubles.
        %
        % The annotation here has been REMAPPED to Allen parcellation_index,
        % the same space data\atlas\annotation_10.nii.gz uses and the only one
        % get_allen_region_mask.m understands. BrainGlobe ships it in Allen
        % structure IDs instead, and the two spaces collide numerically without
        % meaning the same thing -- structure 672 (CP) is index 662. All 686
        % ids translated, no labelled voxel lost. The original is kept beside it
        % as annotation_structureids_original.nii.gz.
        %
        % default_aplims is measured twice per age, by two independent methods
        % that agreed for P20: matching the brain's AP cross-sectional area
        % profile gives [63 559], and regressing the AP centre of mass of 678
        % corresponding regions gives [62 562] (r = 0.998, residual 0.17 mm).
        % Mapping the adult crop across by brain fraction had given [97 566],
        % which is wrong. build_demba_atlas.py runs both and writes the first
        % into aplims.txt, warning if they disagree by more than half a mm.
        %
        % That regression also measures something worth knowing: its slope is
        % 1.797 CCF planes per DeMBA plane, not the 2.000 the voxel sizes imply.
        % DeMBA is about 11% longer in AP than the CCF for the same anatomy.
        % Carey 2025 says why -- the CCFv3 template is rostrocaudally shrunken
        % and the developmental templates are not. A slicethickness tuned
        % against the CCF therefore under-scales against DeMBA by roughly that
        % much. Reproduced by tmp/remap_demba.py and tmp/check_crop.py.

end

% Fail early and clearly rather than deep inside a registration call
if ~exist(atlas.dir, 'dir')
    error('get_atlas: atlas dir not found: %s', atlas.dir);
end
required = {atlas.template_file, atlas.annotation_file};
for k = 1:numel(required)
    fp = fullfile(atlas.dir, required{k});
    if ~exist(fp, 'file')
        error('get_atlas: required atlas file missing: %s', fp);
    end
end

% Every atlas directory holds files with the SAME names, because that is how
% LightSuite finds them. So if two of them are on the MATLAB path at once,
% which() silently picks whichever was added first and a brain can be registered
% to the wrong atlas with nothing in the log to say so. Drop the others here.
% Compare whole path ENTRIES, not substrings. 'atlas_demba_p20' contains
% 'atlas', so a substring test on the path string reports the adult dir as
% present whenever the DeMBA one is, and rmpath then warns about a directory
% that was never there.
% Every built age counts, not just the two that existed when this was written:
% with P16 and P20 both present, a stale entry for one is exactly how a brain
% gets registered to the template of the wrong age.
all_atlas_dirs = {p.atlas};
demba_dirs = dir(fullfile(p.data, 'atlas_demba_p*'));
for k = 1:numel(demba_dirs)
    if demba_dirs(k).isdir
        all_atlas_dirs{end+1} = fullfile(p.data, demba_dirs(k).name); %#ok<AGROW>
    end
end
path_entries = strsplit(path, pathsep);
for k = 1:numel(all_atlas_dirs)
    d = all_atlas_dirs{k};
    if ~strcmpi(d, atlas.dir) && any(strcmpi(path_entries, d))
        rmpath(d);
        fprintf('get_atlas: removed %s from the path so ''%s'' resolves unambiguously.\n', ...
            d, atlas.key);
    end
end
addpath(atlas.dir);

% Cheap proof that the atlas actually in force is the one that was asked for
resolved = which(atlas.template_file);
if ~strcmpi(fileparts(resolved), atlas.dir)
    error(['get_atlas: %s still resolves to\n  %s\ninstead of\n  %s\n' ...
           'Something else put another atlas dir on the path after this call.'], ...
           atlas.template_file, resolved, atlas.dir);
end

end
