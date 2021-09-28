%% Ground Truth

dx = 1;
dy = 1;


x_ = -5e3:dx:5e3; Nx = length(x_);
y_ = -5e3:dy:5e3; Ny = length(y_);

[x, y] = meshgrid(x_, y_);

P_x = x_(1:500); p_y = y_(6000:7000);

Image = meshgrid(zeros(size(x_)), zeros(size(y_)));

Image(1:500, 6000:7000) = 1;


figure; surface(x, y, Image); shading flat; colormap('jet'); colorbar;