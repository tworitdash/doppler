function mat = gauss2d(R, C, sigma, center)
% gsize = size(mat);
% [R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);