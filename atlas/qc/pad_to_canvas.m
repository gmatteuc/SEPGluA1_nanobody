function out = pad_to_canvas(img, canvas_h, canvas_w, fill_value)
%PAD_TO_CANVAS  Centre an image on a canvas filled with one value.
%   out = PAD_TO_CANVAS(img, canvas_h, canvas_w, fill_value) returns a
%   canvas_h x canvas_w array of img's class, filled with fill_value, with img
%   in its centre (the odd pixel of a margin goes to the bottom or the right).
%
%   Run by the atlas QC scripts compare_atlas_regions and
%   compare_atlases_montage.

out = repmat(cast(fill_value, 'like', img), canvas_h, canvas_w);
[h, w] = size(img);
r0 = floor((canvas_h - h) / 2) + 1;
c0 = floor((canvas_w - w) / 2) + 1;
out(r0:r0 + h - 1, c0:c0 + w - 1) = img;

end
