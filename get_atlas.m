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

switch lower(atlas_key)

    case 'ccf'
        atlas.key             = 'ccf';
        atlas.dir             = 'D:\sep_histology\data\atlas';
        atlas.template_file   = 'average_template_10.nii.gz';
        atlas.annotation_file = 'annotation_10.nii.gz';
        atlas.boundary_file   = 'annotation_boundary_10.nii.gz';
        atlas.res_um          = 10;
        atlas.default_aplims  = [180 1079];
        atlas.description     = 'Allen Mouse Brain Common Coordinate Framework v3, 10 um (adult, P56)';

    otherwise
        error('get_atlas: unknown atlas key "%s". Known keys: ''ccf''.', atlas_key);
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

end
