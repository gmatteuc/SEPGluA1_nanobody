% =========================================================================
% Stream sagittal planes (Z = frame) across mice into one 4-D matrix
% 
% =========================================================================

clear; clc;

%% ------------------------------------------------------------------------
% 1) Build the 3-D brain mask (900 × 800 × 1140)
% -------------------------------------------------------------------------
AllenFile = '...\annotation_10.nii.gz'; % Allen CCFv3 NIfTI
AllenVol = niftiread(AllenFile); % (1320×800×1140)
limits = [180 1079]; % crop in AP 
AllenCrop = AllenVol(limits(1):limits(2),:,:); % (900×800×1140)
brainMask = AllenCrop > 0; % logical mask

%% ------------------------------------------------------------------------
% 2) Locate TIFF stacks for each mouse
% -------------------------------------------------------------------------
root = '...\RWS'; %folder containing the group (Naive/RWS/Behavior)
mice = dir(fullfile(root,'MG*')); % folder for each mouse
nMice = numel(mice);
autoFiles = cell(nMice,1);
for m = 1:nMice
    d = fullfile(root, mice(m).name);
    f = dir(fullfile(d,'*_residual.tif*')); % change name depending on the channel we want to load
    assert(~isempty(f), 'No residual TIFF in %s',d); %same 
    autoFiles{m} = fullfile(f(1).folder,f(1).name);
end

%% ------------------------------------------------------------------------
% 3) Inspect first file to get Y, X, Z (Nframes)
% -------------------------------------------------------------------------
info = imfinfo(autoFiles{1});
Nframes = numel(info); % == 1140, sagittal planes
Y = info(1).Height; % 900
X = info(1).Width; % 800
assert(all(size(brainMask)==[Y X Nframes]), 'Mask vs TIFF dimension mismatch');

%% ------------------------------------------------------------------------
% 4) Prepare output directory and disk-backed .mat files
% -------------------------------------------------------------------------
outDir = '...';
if ~exist(outDir,'dir'), mkdir(outDir); end

resMF = matfile(fullfile(outDir,'....mat'),'Writable',true);
resMF.residualVol(Y,X,Nframes,nMice) = single(0); 

%% ------------------------------------------------------------------------
% 5) Stream slices: mask non-brain voxels and write to disk
% -------------------------------------------------------------------------
for m = 1:nMice % Mouse index
    fprintf('Storing mouse %d / %d …\n',m,nMice);
    for f = 1:Nframes % Sagittal slice index (Z)
        sl = single( imread(autoFiles{m}, f) ); % [Y×X] slice
        mask2D = brainMask(:,:,f);
        sl(~mask2D) = NaN; % keep only brain voxels
        % (Y , X , Z , Mouse)
        resMF.residualVol(:,:,f,m) = sl;
    end
end
