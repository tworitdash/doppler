

N = 1024;
PRT = 1e-3;

dr = 10;

[vx, vy, vz, t] = ndgrid(-5e2:dr:5e2, -5e2:dr:5e2, -5e2:dr:5e2, linspace(eps, N*PRT, N)); 

u = 5 * ones(size(vx));
v = zeros(size(vx)); 
w = zeros(size(vx));

px_filter_indx = linspace(1, 25, 25);
py_filter_indx = linspace(1, 25, 25);
pz_filter_indx = linspace(25, 50, 25);

target_x = zeros(size(vx));
target_y = zeros(size(vx));
target_z = zeros(size(vx));


target_x(px_filter_indx) = 1;
target_y(py_filter_indx) = 1;
target_z(py_filter_indx) = 1;

surface(target_z(:, :, 1,1)); shading flat; colormap('jet');