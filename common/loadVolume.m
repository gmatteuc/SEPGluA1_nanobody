function vol = loadVolume(vol_path, ~)
%LOADVOLUME  Read a volume from a multi-page TIFF, or from a folder of TIFFs.
%   vol = LOADVOLUME(vol_path) reads every page of the TIFF file vol_path, or every
%   .tif file of the folder vol_path in natural order (natsortfiles), one plane
%   each, into an H x W x n single array. vol_path may also be a cell holding
%   the path; a second input (the callers pass 1) is ignored.

% the callers pass the path in a cell
if iscell(vol_path)
    vol_path = vol_path{1};
end

if isfolder(vol_path)

    % a folder: one plane per .tif file, in natural order
    files = dir(fullfile(vol_path, '*.tif'));
    files = natsortfiles({files.name});
    n_planes = numel(files);
    info = imfinfo(fullfile(vol_path, files{1}));
    vol = zeros(info.Height, info.Width, n_planes, 'single');
    for z = 1:n_planes
        vol(:, :, z) = single(imread(fullfile(vol_path, files{z})));
    end
else

    % a multi-page TIFF: one plane per page
    info = imfinfo(vol_path);
    n_planes = numel(info);
    vol = zeros(info(1).Height, info(1).Width, n_planes, 'single');
    for z = 1:n_planes
        vol(:, :, z) = single(imread(vol_path, z));
    end
end

end
