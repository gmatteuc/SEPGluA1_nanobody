clear all; clc;
%% Load Data and Atlas

autoPath = '...\chan03_Cy3.tiff'; %autofluorescence registered tiff
nanoPath = '...\chan02_Cy5.tiff'; %nanobody 
allenDir = '...\atlas'; %allen atlas folder including annotation volume and parcellation csv

autoVol = loadVolume({autoPath}, 1); 
nanoVol = loadVolume({nanoPath}, 1);

limits    = [180 1079];  %atlas AP limits
av = niftiread(fullfile(allenDir,'annotation_10.nii.gz'));
av = av(limits(1):limits(2), :, :);   

k         = 2.6;                            % σ‑cut‑off for AF scaling
nbinsHist = 256;


%% Keep the longest run of coronal slices with sufficient signal (drop z with partial slices)
baseMask = av ~= 0;                 % atlas mask 

% Mask with signal in either channel 
epsI = 1e-8;                        
rawValid = baseMask & (autoVol > epsI | nanoVol > epsI);

% Coverage per coronal slice (dim 1): fraction of brain voxels that are valid
numBrainPerSlice = squeeze(sum(sum(baseMask, 2), 3));
numValidPerSlice = squeeze(sum(sum(rawValid, 2), 3));
cov   = numValidPerSlice ./ max(numBrainPerSlice, 1);   % 0..1
covSm = movmean(cov, 5);                                % smooth 
thr   = 0.60;                                           % require ≥60% coverage

% Find the largest contiguous run of slices with cov ≥ thr
good = covSm >= thr;
d = diff([false; good; false]);
runStart = find(d == 1);
runEnd   = find(d == -1) - 1;
if isempty(runStart)
    warning('No slices pass coverage threshold; using all coronal slices.');
    yL = 1; yR = size(baseMask,1);
else
    runLen = runEnd - runStart + 1;
    [~, ix] = max(runLen);         % choose the longest contiguous block
    yL = runStart(ix);
    yR = runEnd(ix);
end

% Final mask
validMask = false(size(baseMask));
validMask(yL:yR,:,:) = rawValid(yL:yR,:,:);

brainMask = validMask;

fprintf('Using largest contiguous coronal block: [%d..%d] (%d slices), thr=%.2f\n', ...
        yL, yR, yR - yL + 1, thr);

%% Use signal within the mask

autoVals=autoVol(brainMask);

nanoVals=nanoVol(brainMask);

%% PERCENTILE‑BASED SCALING 
loPct = 50*(1 - erf(k/sqrt(2)));
hiPct = 100 - loPct;

loValAF = prctile(autoVals, loPct);
hiValAF = prctile(autoVals, hiPct);
autoScaled_vec = (autoVals - loValAF) ./ max(hiValAF - loValAF, eps);
autoScaled_vec = min(max(autoScaled_vec,0),1);

clipFracAF = mean(autoScaled_vec==0 | autoScaled_vec==1)*100;
fprintf('clipped %.3f %% of voxels (k=%.1fσ; lo=%.3f%%, hi=%.3f%%).\n', ...
        clipFracAF, k, loPct, hiPct);

loValNB = prctile(nanoVals, loPct);
hiValNB = prctile(nanoVals, hiPct);

slope     = (hiValAF - loValAF) / max(hiValNB - loValNB, eps);
intercept =  loValAF - slope * loValNB;

nano_inAF      = nanoVals * slope + intercept;  % NB in AF units
nanoScaled_vec = max( (nano_inAF - loValAF) ./ max(hiValAF - loValAF, eps), 0 );  % lower-clip only

nbLowClip = mean(nanoScaled_vec==0)*100;
nbAboveOne = mean(nanoScaled_vec>1)*100;
fprintf('NB scaled: %.3f %% low-clipped; %.3f %% > 1 (by design, no upper clip).\n', ...
        nbLowClip, nbAboveOne);

%% Autofluo removal with Regression using thalamus

parcelinfo = readtable(fullfile(allenDir, ...
    'parcellation_to_parcellation_term_membership.csv'));

roiNames = {
    'paraventricular nucleus of the thalamus',
    'parataenial nucleus',
    'nucleus of reuniens',
    'rhomboid nucleus',
    'xiphoid thalamic nucleus',
    'central medial nucleus',
    'paracentral nucleus',
    'central lateral nucleus',
    % 'ventral anterior-lateral complex of the thalamus',   
    % 'submedial nucleus'
     };

isStruct  = strcmp(parcelinfo.parcellation_term_set_name,'structure');
namesLow  = lower(parcelinfo.parcellation_term_name);
keep_ids  = parcelinfo.parcellation_index( isStruct & ismember(namesLow, roiNames) );

roiMask = ismember(av, keep_ids);
regMask = roiMask & brainMask;        
regMask2d=regMask(brainMask);

autoScaled_vec=double(autoScaled_vec); %cast to double for the regression
nanoScaled_vec=double(nanoScaled_vec);

% predictor & response inside regMask
X = autoScaled_vec(regMask2d);          % AF 0–1
Y = nanoScaled_vec(regMask2d);        % NB rescaled

b = robustfit(X, Y, 'ols');             
fprintf('New slope = %.3f, intercept = %.3f\n', b(2), b(1));

% voxel-wise estimate and residual b
bgEstFull  = b(1) + b(2) * autoScaled_vec;
residual = nanoScaled_vec - bgEstFull;      

%% Scatter and VPlot

X = autoScaled_vec(regMask2d);   X = X(:);   % AF 0–1
Y = nanoScaled_vec(regMask2d); Y = Y(:);   % NB 


% Visual check on random 5 % subsample

idx      = randperm(numel(X), round(0.05*numel(X)));
Xsub     = X(idx);
Ysub     = Y(idx);
Rsub        = residual(regMask2d);  Rsub = Rsub(idx);% residual

figure;
scatter(Xsub, Ysub, 20, 'filled','MarkerFaceAlpha',0.2); hold on
xLine    = linspace(0, 1, 200);
plot(xLine, b(1)+b(2)*xLine, 'k-', 'LineWidth',2);
xlabel('AF scaled (0–1)'); ylabel('NB scaled (0–1)');
title(sprintf('ROI regression: slope %.3f  intercept %.3f', b(2), b(1)));
axis square; xlim([0 1]); ylim([0 1]);


figure;
    violinplot([Xsub, Ysub, Rsub], ...
               {'AF (X)', 'NB (y)', 'Residual'});
    ylabel('Scaled intensity');
    title(sprintf(['Intensity distributions in %d ROI voxels ' ...
                   '(violin, k = %.1f\\sigma)'], numel(idx), k));

%% WRITE TIFF

% 1) Rebuild a 3-D float volume with NaNs outside the mask
vol3d = nan(size(brainMask), 'single');
vol3d(brainMask) = single(residual);   

% 2) Write as BigTIFF (32-bit float, LZW compressed)
outFile = '';
[nX,nY,nZ] = size(vol3d);

t = Tiff(outFile, 'w8');   % "w8" = BigTIFF
tag.ImageLength         = nX;                                   % rows  = dim-1 (coronal)
tag.ImageWidth          = nY;                                   % cols  = dim-2
tag.SampleFormat        = Tiff.SampleFormat.IEEEFP;             % float32
tag.BitsPerSample       = 32;
tag.SamplesPerPixel     = 1;
tag.Photometric         = Tiff.Photometric.MinIsBlack;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression         = Tiff.Compression.LZW;
tag.RowsPerStrip        = nX;                                   

for k = 1:nZ
    setTag(t, tag)
    slice = vol3d(:,:,k);              % keep orientation as in memory
    t.write(slice)
    if k < nZ, t.writeDirectory(); end
end
t.close();

fprintf('Wrote BigTIFF: %s  (%dx%dx%d float)\n', outFile, nX, nY, nZ);

%also save the mask
maskFile = strrep(outFile,'.tif','_validMask.tif');
tm = Tiff(maskFile, 'w8');
tag.BitsPerSample = 8; tag.SampleFormat = Tiff.SampleFormat.UInt;
for k = 1:nZ
    setTag(tm, tag)
    tm.write(uint8(brainMask(:,:,k))*255);
    if k < nZ, tm.writeDirectory(); end
end
tm.close();
