function [row, col] = lr_atlas_boundaries(atlas_vol, j)
%LR_ATLAS_BOUNDARIES  Region boundaries of one AP plane of the atlas, for the videos.
%   [row, col] = LR_ATLAS_BOUNDARIES(atlas_vol, j) returns the DV and ML
%   positions on plane j of atlas_vol (AP x DV x ML) where the annotation
%   changes along ML, inside regions with an id above 1.
%
%   Run by write_lr_video, write_lr_video_surpmask and
%   write_lr_video_surpmask_rolling.

atlasim = squeeze(atlas_vol(j, :, :));
atlasim = single(atlasim);
av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
[row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

end
