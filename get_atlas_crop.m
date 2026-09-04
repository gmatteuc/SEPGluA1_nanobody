function A = get_atlas_crop(atlas_key)
%GET_ATLAS_CROP The annotation volume on the grid of the registered volumes.
%
%   A = get_atlas_crop('ccf')         adults
%   A = get_atlas_crop('demba_p20')   the P20 cohort
%
% Returns a struct:
%   annot      cropped annotation volume, AP x DV x ML, on the same grid as
%              volume_registered\chan0X_*.tiff for that cohort
%   brainMask  annot > 0
%   n_ap       number of AP planes (900 for the adults, 994 for P20)
%   ap_scale   n_ap / 900: what an adult AP index has to be multiplied by to
%              land at the same relative position in this cohort's volume
%   key, res_um, aplims, csv_dir
%
% Why this exists. P5 to P10 all opened annotation_10.nii.gz and cropped it
% to [180 1079] inline, which is the adult atlas and nothing else. The young
% brains are registered to DeMBA, and their registered volumes come out on a
% different grid: LightSuite writes them at twice the 20 um registration
% grid, so the adults land on the CCF 10 um grid (900 x 800 x 1140) and the
% P20 brains on a 10 um-equivalent grid of DeMBA space (994 x 800 x 1140,
% i.e. the 20 um crop [63 559] upsampled by two). DV and ML are the same size
% in both, only the AP length differs. This helper produces the matching
% annotation for either, so the analysis scripts can stop knowing which.
%
% For 'ccf' the returned annot is exactly AllenVol(180:1079, :, :) -- the
% same array, to the byte, that the scripts used to build themselves. That
% is what keeps the adult results identical.
%
% The label ids in the DeMBA annotation were remapped to Allen
% parcellation_index when the atlas was built, so the ontology CSVs in the
% CCF folder (csv_dir) apply to both; get_allen_region_mask takes csv_dir.

paths = get_paths();
atlas = get_atlas(atlas_key);

vol = niftiread(fullfile(atlas.dir, atlas.annotation_file));
lims = atlas.default_aplims;
annot = vol(lims(1):lims(2), :, :);
clear vol

% The registered volumes sit on a 10 um grid whatever the atlas resolution,
% so an atlas coarser than 10 um is brought up to it. Nearest neighbour, so
% labels stay labels.
factor = atlas.res_um / 10;
if factor ~= 1
    annot = imresize3(annot, factor, 'Method', 'nearest');
end

A = struct();
A.key       = atlas.key;
A.res_um    = atlas.res_um;
A.aplims    = lims;
A.csv_dir   = paths.atlas;                 % ontology CSVs live with the CCF
A.annot     = annot;
A.brainMask = annot > 0;
A.n_ap      = size(annot, 1);
A.ap_scale  = A.n_ap / 900;                % adults: 900 planes, scale 1
end
