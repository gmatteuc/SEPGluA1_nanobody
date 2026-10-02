function m = region_overlap_matrix(csv_dir, av_ref, av_test, region_names)
%REGION_OVERLAP_MATRIX  How much of each named region one annotation gives each other.
%   m = REGION_OVERLAP_MATRIX(csv_dir, av_ref, av_test, region_names) returns
%   the n x n matrix, n = numel(region_names), whose element (i, j) is the
%   fraction of region i of the annotation av_ref (rows) that the annotation
%   av_test calls region j (columns); a row stays 0 when av_ref has none of
%   that region. A region is the named Allen region with its descendants
%   (get_allen_region_mask, with the ontology files in csv_dir), within the
%   voxels its annotation labels.
%
%   Run by the atlas QC scripts check_demba_to_allen and compare_atlas_regions.

n = numel(region_names);
mask_ref  = cell(1, n);
mask_test = cell(1, n);
for i = 1:n
    mask_ref{i}  = get_allen_region_mask(csv_dir, av_ref,  region_names(i), av_ref > 0);
    mask_test{i} = get_allen_region_mask(csv_dir, av_test, region_names(i), av_test > 0);
end
m = zeros(n);
for i = 1:n
    denom = nnz(mask_ref{i});
    if denom == 0
        continue
    end
    for j = 1:n
        m(i, j) = nnz(mask_ref{i} & mask_test{j}) / denom;
    end
end

end
