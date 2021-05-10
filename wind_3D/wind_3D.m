%% Wind field in 3D space
nx = 1000; 
ny = 1000;
nz = 1000;


vx = linspace(eps, 5e3, nx); % positions of x, y and z in space
vy = linspace(eps, 5e3, ny);
vz = linspace(eps, 5e3, nz);

[x, y, z] = ndgrid(vx, vy, vz); % ndgrid of position vectors

ux = 5;

u = ux  .* ones(size(x));
v = zeros(size(x));
w = zeros(size(x));

%%

quiver3(x, y, z, u, v, w);

