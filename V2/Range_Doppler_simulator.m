

N = 1024;
PRT = 1e-3;

[vx, vy, vz, t] = ndgrid(linspace(-5e3, 5e3, 100), linspace(-5e3, 5e3, 100), linspace(-5e3, 5e3, 100), linspace(eps, N*PRT, N)); 

u = 5 * ones(size(x));
v = zeros(x); 
w = zeros(x);

% [px, py, pz, t] = 