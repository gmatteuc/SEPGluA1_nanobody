function R = sep_compare_figure_images(f_a, f_b, P)
%SEP_COMPARE_FIGURE_IMAGES  Whether a regenerated figure is the same as the original.
%   R = SEP_COMPARE_FIGURE_IMAGES(f_a, f_b, P) compares two rendered figures
%   (PNG or JPG) and returns a verdict in R.verdict and the numbers behind it.
%
%   Two renderings of the same figure can differ where anti-aliasing puts an
%   edge or a glyph a fraction of a pixel elsewhere. Anti-aliasing colours a
%   pixel with a mix of the colours drawn over it, so moving the drawing by
%   less than a pixel changes the mix, never the colours: the new pixel is a
%   weighted average of the original's colours within 1 px of it. A changed
%   pixel therefore counts as rendering only if, within noise_tol in every
%   colour channel, it is a mix of at most three of the original's colours
%   in its 3x3 neighbourhood (two for an edge, three where edges meet). Any
%   other changed pixel is a real change and makes the figure DIFFERENT: a
%   colour changed at the same brightness, a changed digit, ink where there
%   was none. The colour channels are compared separately, never as grey.
%   Images of different sizes are DIFFERENT, never resized.
%
%   noise_tol allows for rounding, and for a line that meets another
%   element: there a pixel can show the line's own colour, which no pixel of
%   the original around it shows at full strength (6 levels off in the
%   tests).
%
%   What this cannot see: next to an edge, a pixel may take any mix of the
%   colours on either side, so a change confined to the anti-aliased rim of
%   ink passes. Changed text, lines and areas reach pixels whose
%   neighbourhood held one colour only, but a small change reaches few: a 9
%   turned into an 8 at font size 10 left 2 such pixels in the tests. JPEG
%   compression error is no mix: a JPEG drawn a fraction of a pixel
%   elsewhere comes out DIFFERENT, so compare figures saved as PNG.
%
%   Options (P):
%     noise_tol   distance from a mix still counted as
%                 rendering, per channel on 0-255                  (8)
%     roi         [r0 r1] fraction of height judged, for figures
%                 whose titles differ by design                    ([0 1])
%     diff_image  path of a difference image, real changes red
%                 on the original                                  ('')
%
%   R: ok, verdict, n_changed (pixels that differ at all), n_real (changed
%   pixels that no mix explains), max_dist (largest distance of a changed
%   pixel from a mix, 0-255: how close a match came to noise_tol), size_a,
%   size_b (rows, columns, colour channels).
%
%   See also SEP_COMPARE_OUTPUTS.

if nargin < 3
    P = struct();
end

% defaults
d.noise_tol  = 8;
d.roi        = [0 1];
d.diff_image = '';

% fill missing or empty options with the defaults
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(P,f{i}) || isempty(P.(f{i}))
        P.(f{i}) = d.(f{i});
    end
end

R = struct('ok',false, 'verdict','MISSING', 'n_changed',NaN, 'n_real',NaN, ...
           'max_dist',NaN, 'size_a','', 'size_b','');

% stop if either file is missing
if exist(f_a,'file') ~= 2
    R.verdict = 'ORIGINAL MISSING';
    return;
end
if exist(f_b,'file') ~= 2
    R.verdict = 'REGENERATED MISSING';
    return;
end

a = read_rgb(f_a);
b = read_rgb(f_b);
R.size_a = mat2str(size(a));
R.size_b = mat2str(size(b));

% a different pixel size is a different figure: resizing would blur it into
% a match
if ~isequal(size(a), size(b))
    R.verdict = 'DIFFERENT';
    return;
end

% optional band restriction, for figures whose titles differ by design
[n_rows, n_cols, n_channels] = size(a);
band = false(n_rows, n_cols);
r0 = max(1, round(P.roi(1)*n_rows) + 1);
r1 = min(n_rows, round(P.roi(2)*n_rows));
band(r0:r1, :) = true;
if ~any(band(:))
    R.verdict = 'EMPTY ROI';
    return;
end

% the changed pixels, one per row
changed = find(band & any(a ~= b, 3));
a_list = reshape(a, [], n_channels);
b_list = reshape(b, [], n_channels);

% how far each changed pixel is from the closest mix of the original's
% colours around it; in chunks, so that an image that changed everywhere
% never needs nine copies of itself at once
D = zeros(numel(changed), 1);
chunk = 2^18;
for i0 = 1:chunk:numel(changed)
    part = i0:min(numel(changed), i0 + chunk - 1);
    around = neighbour_colours(a_list, changed(part), n_rows, n_cols);
    D(part) = mix_distance(around, b_list(changed(part), :));
end
unexplained = changed(D > P.noise_tol);

R.n_changed = numel(changed);
R.n_real    = numel(unexplained);
R.max_dist  = max([0; D]);
R.ok        = R.n_real == 0;
if R.ok
    R.verdict = 'MATCH';
else
    R.verdict = 'DIFFERENT';
end

% difference image: real changes in red on the original, grown by a pixel
% so that a single one can be seen
if ~isempty(P.diff_image)
    M = false(n_rows, n_cols);
    M(unexplained) = true;
    M = conv2(double(M), ones(3), 'same') > 0;
    grey = rgb2gray(uint8(a));
    r = grey;
    g = grey;
    bl = grey;
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

function I = read_rgb(f)
% An image as RGB on 0-255, whatever its storage: grey, 8 or 16 bit, or
% indexed (through its colour map, so a changed map is a changed image).

[I, map] = imread(f);
if ~isempty(map)
    I = ind2rgb(I, map);
end
I = 255 * im2double(I);
if size(I,3) == 1
    I = repmat(I, 1, 1, 3);
end
end

function colours = neighbour_colours(a_list, pixels, n_rows, n_cols)
% The original's colours in the 3x3 neighbourhood of each pixel (linear
% indices): one cell per neighbour, the pixel itself included, one row per
% pixel. The border is repeated at the edges of the image.

[rr, cc] = ind2sub([n_rows n_cols], pixels);
colours = cell(1, 9);
k = 0;
for dy = -1:1
    for dx = -1:1
        k = k + 1;
        rows = min(max(rr + dy, 1), n_rows);
        cols = min(max(cc + dx, 1), n_cols);
        colours{k} = a_list(sub2ind([n_rows n_cols], rows, cols), :);
    end
end
end

function D = mix_distance(colours, x)
% For each row of x, how far that colour is from the closest mix of at most
% three of the colours in the same row of colours{1}, colours{2}, ...: the
% largest channel difference. The mixes of three colours fill a triangle,
% whose sides are the mixes of two; a point outside a triangle is nearest to
% one of its sides, so the triangles are only needed for their inside.

n = numel(colours);
D = inf(size(x, 1), 1);

% mixes of two colours: the segment from one to the other
pairs = nchoosek(1:n, 2);
for k = 1:size(pairs, 1)
    c0 = colours{pairs(k, 1)};
    u = colours{pairs(k, 2)} - c0;
    v = x - c0;
    uu = sum(u.^2, 2);

    % the closest point of the segment, by projection; where both colours
    % are the same, the only mix is that colour
    t = sum(u .* v, 2) ./ uu;
    t(uu == 0) = 0;
    t = min(max(t, 0), 1);
    D = min(D, max(abs(v - t .* u), [], 2));
end

% mixes of three colours: inside the triangle, c0 + s*u1 + t*u2 with s and
% t at least 0 and s + t at most 1 (least squares, by the normal equations)
triples = nchoosek(1:n, 3);
for k = 1:size(triples, 1)
    c0 = colours{triples(k, 1)};
    u1 = colours{triples(k, 2)} - c0;
    u2 = colours{triples(k, 3)} - c0;
    v = x - c0;
    g11 = sum(u1.^2, 2);
    g12 = sum(u1 .* u2, 2);
    g22 = sum(u2.^2, 2);
    h1 = sum(u1 .* v, 2);
    h2 = sum(u2 .* v, 2);
    denom = g11 .* g22 - g12.^2;
    s = (g22 .* h1 - g12 .* h2) ./ denom;
    t = (g11 .* h2 - g12 .* h1) ./ denom;

    % three colours on one line make no triangle: its sides cover them
    inside = denom > 0 & s >= 0 & t >= 0 & s + t <= 1;
    d = max(abs(v - s .* u1 - t .* u2), [], 2);
    D(inside) = min(D(inside), d(inside));
end
end
