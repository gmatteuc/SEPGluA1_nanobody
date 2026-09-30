function R = sep_compare_figure_images(f_a, f_b, P)
%SEP_COMPARE_FIGURE_IMAGES  Whether a regenerated figure is the same as the original.
%   R = SEP_COMPARE_FIGURE_IMAGES(f_a, f_b, P) compares two rendered figures
%   (PNG or JPG) and returns a verdict in R.verdict and the numbers behind it.
%
%   Counting differing pixels never resolves: two MATLAB sessions rasterise
%   glyphs, anti-aliasing and JPEG ringing differently, 2-3% of pixels. Those
%   differences sit on edges and move ink by a pixel or two; a wrong number
%   changes the flat interior of something (a cell, a map, a band). So only the
%   flat part is judged, with the same slack:
%     edge  within slack_px of a line or glyph; reported, never judged
%     flat  everything else; this is the verdict
%   Validated on the imaging repository's figures: a different figure of the
%   same kind scores 4-8% flat, a true match 0.0002-0.0004%, so the test can
%   fail.
%
%   Options (P):
%     diff_tol    pixel difference that counts as visible, 0-255   (12)
%     slack_px    how far ink may move before it is content        (2)
%     edge_grad   gradient above which a pixel is line or glyph    (8)
%     flat_pass   flat % allowed before DIFFERENT                  (0.01)
%     roi         [r0 r1] fraction of height judged, for figures
%                 whose titles differ by design                    ([0 1])
%     diff_image  path of a difference image, flat differences red
%                 on the original                                  ('')
%
%   R: ok, verdict, flat_pct, flat_max, flat_px, edge_pct, size_a, size_b,
%   resized.
%
%   See also SEP_COMPARE_OUTPUTS.

if nargin < 3
    P = struct();
end

% defaults
d.diff_tol   = 12;
d.slack_px   = 2;
d.edge_grad  = 8;
d.flat_pass  = 0.01;
d.roi        = [0 1];
d.diff_image = '';

% fill missing or empty options with the defaults
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(P,f{i}) || isempty(P.(f{i}))
        P.(f{i}) = d.(f{i});
    end
end

R = struct('ok',false, 'verdict','MISSING', 'flat_pct',NaN, 'flat_max',NaN, ...
           'flat_px',0, 'edge_pct',NaN, 'size_a','', 'size_b','', 'resized',false);

% stop if either file is missing
if exist(f_a,'file') ~= 2
    R.verdict = 'ORIGINAL MISSING';
    return;
end
if exist(f_b,'file') ~= 2
    R.verdict = 'REGENERATED MISSING';
    return;
end

A = imread(f_a);
B = imread(f_b);
R.size_a = mat2str(size(A,1:2));
R.size_b = mat2str(size(B,1:2));

% a different pixel size means a different window: resize, but record it
if ~isequal(size(A,1:2), size(B,1:2))
    B = imresize(B, size(A,1:2));
    R.resized = true;
end

% compare in grey levels
a = double(gray_safe(A));
b = double(gray_safe(B));

% smallest difference achievable by letting the ink slide a little
D = slack_diff(a, b, P.slack_px);

% edges of the original, dilated by the same slack: lines, glyphs and the
% anti-aliased skirt around them
[gy, gx] = gradient(a);
edge = dilate(hypot(gy,gx) > P.edge_grad, P.slack_px);

% optional band restriction, for figures whose titles differ by design
band = false(size(a));
r0 = max(1, round(P.roi(1)*size(a,1)) + 1);
r1 = min(size(a,1), round(P.roi(2)*size(a,1)));
band(r0:r1, :) = true;

% judge the flat part only
flat = ~edge & band;
if ~any(flat(:))
    R.verdict = 'NO FLAT AREA';
    return;
end

R.flat_px  = nnz(flat);
R.flat_pct = 100 * mean(D(flat) > P.diff_tol);
R.flat_max = max(D(flat));
R.edge_pct = 100 * mean(D(edge & band) > P.diff_tol);
R.ok       = R.flat_pct <= P.flat_pass;
if R.ok
    R.verdict = 'MATCH';
else
    R.verdict = 'DIFFERENT';
end

% difference image: flat differences in red on the original
if ~isempty(P.diff_image)
    M   = dilate(flat & (D > P.diff_tol), 1);
    rgb = repmat(uint8(a), 1, 1, 3);
    r = rgb(:,:,1);
    g = rgb(:,:,2);
    bl = rgb(:,:,3);
    r(M) = 255;
    g(M) = 0;
    bl(M) = 0;

    % save, creating the folder if needed
    dd = fileparts(P.diff_image);
    if ~isempty(dd) && ~exist(dd,'dir')
        mkdir(dd);
    end
    imwrite(cat(3,r,g,bl), P.diff_image);
end
end

% ===== Local functions =====

function D = slack_diff(a, b, r)
% Smallest |a-b| over shifts of b within +-r px (circular shifts only wrap the
% blank margin).

D = inf(size(a));
for dy = -r:r
    for dx = -r:r
        D = min(D, abs(a - circshift(b, [dy dx])));
    end
end
end

function M = dilate(M, r)
% Binary dilation by a (2r+1) square.

out = M;
for dy = -r:r
    for dx = -r:r
        out = out | circshift(M, [dy dx]);
    end
end
M = out;
end

function g = gray_safe(I)
% Greyscale of an RGB or grey image.

if ndims(I) == 3
    g = rgb2gray(I);
else
    g = I;
end
end
