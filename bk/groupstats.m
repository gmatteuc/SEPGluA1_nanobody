clear all; clc;

% List of tiffs
files = {
  '...\....tiff'
  '...\....tiff'
  '...\....tiff'
   '...\....tiff'
    '...\....tiff'
};

% Get stack size from the first file
info = imfinfo(files{1});
nZ = numel(info); nX = info(1).Height; nY = info(1).Width;

sumVol   = zeros(nX,nY,nZ,'double');
sumSqVol = zeros(nX,nY,nZ,'double');
cntVol   = zeros(nX,nY,nZ,'double');


for f = 1:numel(files)
    t = Tiff(files{f}, 'r');
    try
        for k = 1:nZ
            setDirectory(t, k);
            sl  = single(t.read());          % slice (float), NaN outside valid
            fin = isfinite(sl);              % valid coverage for this mouse

          
            sumVol(:,:,k)   = sumVol(:,:,k)   + sl .* single(fin);
            sumSqVol(:,:,k) = sumSqVol(:,:,k) + (sl.^2) .* single(fin);
            cntVol(:,:,k)   = cntVol(:,:,k)   + uint8(fin);
        end
    catch ME
        t.close();
        rethrow(ME);
    end
    t.close();
end


% mean / SD / SEM 
meanVol = sumVol ./ max(single(cntVol), 1);
meanVol(cntVol==0) = NaN;

varVol = (sumSqVol - single(cntVol).*meanVol.^2) ./ max(single(cntVol)-1, 1);
varVol(cntVol<=1) = NaN;                
stdVol = sqrt(varVol);
semVol = stdVol ./ sqrt(max(single(cntVol),1));

outMat = fullfile('...', '...mat');
save(outMat, 'meanVol','varVol','stdVol','semVol','-v7.3');
fprintf('Saved MAT: %s\n', outMat);

clear all; clc;