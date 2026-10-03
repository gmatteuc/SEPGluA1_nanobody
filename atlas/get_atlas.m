function atlas = get_atlas(atlas_key)
%GET_ATLAS  Location and parameters of a reference atlas, by key.
%   atlas = GET_ATLAS() returns the default, 'ccf'.
%   atlas = GET_ATLAS('ccf') returns the Allen Mouse Brain CCFv3, 10 um.
%   atlas = GET_ATLAS('demba_p20') returns DeMBA at P20, 20 um; any age N
%   built by build_demba_atlas.py is 'demba_pN'.
%
%   It also puts that atlas's folder on the MATLAB path and takes every other
%   atlas folder off it: every atlas folder holds files of the same names,
%   because that is how LightSuite finds them (which()), so with two of them
%   on the path a brain could be registered to the wrong atlas, or the wrong
%   age, with nothing in the log to say so. It stops if a required file is
%   missing, rather than deep inside a registration call.
%
%   atlas is a struct with the fields:
%     key              short identifier used in output tags
%     dir              folder holding the atlas volumes
%     template_file    grayscale average template (used for registration)
%     annotation_file  region-label volume (used for ROI masks)
%     boundary_file    region-boundary volume (used for overlays; '' for DeMBA)
%     res_um           isotropic voxel size in micrometres
%     age_days         postnatal age this atlas represents (56 = adult). An
%                      age-matched atlas is valid only for mice of that age:
%                      a P36 brain does not belong on the P20 template.
%     default_aplims   AP crop of this atlas, in its own planes
%     description      a note for people
%
%   'ccf' gives exactly the folder and files every existing result was
%   produced with, the atlas folder under the data root.
%
%   default_aplims: the registration reads the AP crop from each mouse's
%   local_settings.txt (atlasaplims), not from here, and the drivers from
%   run_extract_and_center to run_register_to_atlas must keep using
%   sliceinfo.atlasaplims so existing results stay as they are. get_atlas_crop
%   reads default_aplims.
%
%   The DeMBA entries are built to one recipe by build_demba_atlas.py and
%   differ only in the age and in the AP crop measured for it. Three things
%   about them are deliberate and easy to get wrong:
%   - The files are named *_10.nii.gz but hold 20 um data. LightSuite finds
%     the atlas with which('average_template_10.nii.gz') in fourteen
%     different files, so the name is forced by the vendored code. The real
%     resolution is res_um here and px_atlas in local_settings.txt, and both
%     must say 20 for a young brain or the AP scale silently doubles.
%   - The annotation is remapped to Allen parcellation_index, the space of
%     data\atlas\annotation_10.nii.gz and the only one get_allen_region_mask
%     understands. BrainGlobe ships it in Allen structure ids, and the two
%     spaces collide numerically without meaning the same thing: structure
%     672 (CP) is index 662. All 686 ids are translated, no labelled voxel is
%     lost, and the original is kept beside it as
%     annotation_structureids_original.nii.gz.
%   - default_aplims comes from aplims.txt beside the volumes, measured when
%     the folder is built, so it cannot drift away from its atlas. Two
%     independent methods agreed for P20: matching the brain's AP
%     cross-sectional area profile gives [63 559], and regressing the AP
%     centre of mass of 678 corresponding regions gives [62 562] (r = 0.998,
%     residual 0.17 mm). Mapping the adult crop across by brain fraction gives
%     [97 566], which is wrong. build_demba_atlas.py runs both and writes the
%     first into aplims.txt, warning if they disagree by more than half a mm.
%
%   That regression also measures the AP stretch: its slope is 1.797 CCF
%   planes per DeMBA plane, not the 2.000 the voxel sizes imply, so DeMBA is
%   about 11% longer in AP than the CCF for the same anatomy. The CCFv3
%   template is rostrocaudally shrunken and the developmental templates are
%   not (Carey 2025). A slicethickness tuned against the CCF therefore
%   under-scales against DeMBA by roughly that much. build_demba_atlas.py
%   measures these numbers for each age it builds (the remap, both crops and
%   the slope, in the folder's source.txt); atlas\qc\atlas_diagnostics
%   measures them again.
%
%   Another atlas gets a key of its own here, pointing at its template and
%   annotation volumes. A comparison across atlases needs both to resolve
%   into a common space (DevCCF ships CCF-linked labels for this); never swap
%   atlases for a comparison without that mapping.

%% Resolve the key

if nargin < 1 || isempty(atlas_key)
    atlas_key = 'ccf';
end

paths = get_paths();

switch lower(atlas_key)

    case 'ccf'
        atlas.key             = 'ccf';
        atlas.dir             = paths.atlas;
        atlas.template_file   = 'average_template_10.nii.gz';
        atlas.annotation_file = 'annotation_10.nii.gz';
        atlas.boundary_file   = 'annotation_boundary_10.nii.gz';
        atlas.res_um          = 10;
        atlas.age_days        = 56;
        atlas.default_aplims  = [180 1079];
        atlas.description     = ...
            'Allen Mouse Brain Common Coordinate Framework v3, 10 um (adult, P56)';

    otherwise

        % a DeMBA age, built by build_demba_atlas.py
        atlas = demba_atlas(atlas_key, paths);

end

%% Check the files

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

%% Put this atlas alone on the path

% this atlas on the path, the others off it
put_alone_on_path(paths, atlas);

end

% ===== Local functions =====

function atlas = demba_atlas(atlas_key, paths)
% The entry of a DeMBA age, its AP crop read from aplims.txt; stops if the age
% is not built.

% any DeMBA age ('demba_p16', 'demba_p20', ...): one entry for all, since
% the folders differ only in the age and its crop (see the help)
tok = regexp(lower(atlas_key), '^demba_p(\d+)$', 'tokens', 'once');
if isempty(tok)
    error(['get_atlas: unknown atlas key "%s". Known keys: ''ccf'' and ' ...
           '''demba_pN'' for any age N built by build_demba_atlas.py.'], atlas_key);
end
age = str2double(tok{1});
atlas.key             = sprintf('demba_p%d', age);
atlas.dir             = fullfile(paths.data, sprintf('atlas_demba_p%d', age));
atlas.template_file   = 'average_template_10.nii.gz';
atlas.annotation_file = 'annotation_10.nii.gz';
atlas.boundary_file   = '';
atlas.res_um          = 20;
atlas.age_days        = age;
atlas.description     = sprintf(['DeMBA P%d (Carey 2025), Allen CCFv3 labels, ' ...
                                 '20 um isotropic'], age);

% an age whose folder has not been built is an error, not a fall-back
% onto a neighbouring age
if ~exist(atlas.dir, 'dir')
    error(['get_atlas: no atlas built for P%d.\n  %s does not exist.\n' ...
           'Build it first:  tools\\venv_atlas\\Scripts\\python.exe atlas\\build_demba_atlas.py %d'], ...
           age, atlas.dir, age);
end

% the AP crop, measured when the folder was built and stored beside the
% volumes (the two methods are in build_demba_atlas.py)
aplims_file = fullfile(atlas.dir, 'aplims.txt');
if ~exist(aplims_file, 'file')
    error(['get_atlas: %s is missing. Re-run atlas\\build_demba_atlas.py %d, or write the ' ...
           'two AP crop planes into that file.'], aplims_file, age);
end
lims = sscanf(fileread(aplims_file), '%d')';
if numel(lims) ~= 2 || lims(2) <= lims(1)
    error('get_atlas: %s should hold two increasing plane numbers, found "%s".', ...
        aplims_file, strtrim(fileread(aplims_file)));
end
atlas.default_aplims  = lims;

end

function put_alone_on_path(paths, atlas)
% The atlas folder on the path and every other one off it, checked by which.

% every atlas folder: the CCF one and every DeMBA age built
all_atlas_dirs = {paths.atlas};
demba_dirs = dir(fullfile(paths.data, 'atlas_demba_p*'));
for k = 1:numel(demba_dirs)
    if demba_dirs(k).isdir
        all_atlas_dirs{end+1} = fullfile(paths.data, demba_dirs(k).name); %#ok<AGROW>
    end
end

% take the others off the path, comparing whole path entries rather than
% substrings ('atlas_demba_p20' contains 'atlas')
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

% check that the atlas in force is the one asked for
resolved = which(atlas.template_file);
if ~strcmpi(fileparts(resolved), atlas.dir)
    error(['get_atlas: %s still resolves to\n  %s\ninstead of\n  %s\n' ...
           'Something else put another atlas dir on the path after this call.'], ...
           atlas.template_file, resolved, atlas.dir);
end

end
