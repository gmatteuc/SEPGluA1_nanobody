function c = get_color2color_colormap(color1, color2)
%GET_COLOR2COLOR_COLORMAP  A diverging colormap from one colour through white to another.
%   c = GET_COLOR2COLOR_COLORMAP(color1, color2) returns a 256 x 3 colormap
%   that runs linearly from the RGB triplet color1 to white over its first
%   half, and from white to color2 over its second. The difference maps use
%   GET_COLOR2COLOR_COLORMAP([0 0 1], [1 0 0]), blue to red, with symmetric
%   colour limits so that zero is white.

% intermediate colour
intermediate_color = [1, 1, 1];

% number of steps of the colormap
if size(gray, 1) ~= 256
    m = 256;
else
    m = size(gray, 1);
end

% the colormap in two halves
m1 = m ./ 2;

% red from color1 to the intermediate colour to color2
r = [linspace(color1(1), intermediate_color(1), m1), ...
    linspace(intermediate_color(1), color2(1), m1)];

% green, the same way
g = [linspace(color1(2), intermediate_color(2), m1), ...
    linspace(intermediate_color(2), color2(2), m1)];

% blue, the same way
b = [linspace(color1(3), intermediate_color(3), m1), ...
    linspace(intermediate_color(3), color2(3), m1)];

% one row per step
c = [r; g; b]';

end
