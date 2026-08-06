function c = get_color2color_colormap(color1,color2)

% Set intermediate color
intermediate_color=[1,1,1];
% Set number of steps of colormap
if size(gray,1)~=256
    m = 256;
else
    m = size(gray,1);
end
% Draw colormap in two steps
m1 = m./2;
% From red of color1 to red of intermediate color to red of color2
r = [linspace(color1(1),intermediate_color(1),m1),linspace(intermediate_color(1),color2(1),m1)];
% From green of color1 to green of intermediate color to green of color2
g = [linspace(color1(2),intermediate_color(2),m1),linspace(intermediate_color(2),color2(2),m1)];
% From blue of color1 to blue of intermediate color to blue of color2
b = [linspace(color1(3),intermediate_color(3),m1),linspace(intermediate_color(3),color2(3),m1)];
% Assign values to colormap
c = [r;g;b]';