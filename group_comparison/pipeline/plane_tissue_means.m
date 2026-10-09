function med_data_4d = plane_tissue_means(data_4d, bkgmask_4d)
%PLANE_TISSUE_MEANS  Each mouse's mean tissue intensity per plane.
%   med_data_4d = PLANE_TISSUE_MEANS(data_4d, bkgmask_4d) returns a planes x
%   mice array: the mean of each plane of each mouse outside its background
%   mask (true for background), NaN left out, and NaN for a plane with no
%   tissue. data_4d and bkgmask_4d are AP x DV x ML x mouse; a 3D volume is
%   one mouse. These are the profiles align_exp_to_ctrl fits its line on.
%   Used by group_differences and per_mouse_region_values.

% one value per plane and mouse
med_data_4d = nan(size(data_4d, 1), size(data_4d, 4));
total_slices = size(data_4d, 1);
for iii = 1:size(data_4d, 4)
    fprintf('  Processing Mouse %d ...\n', iii);
    for slice_idx_loop = 1:total_slices

        % the mean outside the background mask, without the NaN of voxels no
        % section reached; NaN for a plane with no tissue
        img_data = squeeze(data_4d(slice_idx_loop, :, :, iii));
        bg_mask = squeeze(bkgmask_4d(slice_idx_loop, :, :, iii));
        med_data_4d(slice_idx_loop, iii) = nanmean(img_data(~bg_mask)); %#ok<NANMEAN>
    end
end
end
