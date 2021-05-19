%% Wind field in 3D space
clear;
close all;

nx = 1000; 
ny = 1000;
nz = 1000;


vx = linspace(-5e2, 5e2, nx); % positions of x, y and z in space
vy = linspace(-5e2, 5e2, ny);
vz = linspace(-5e2, 5e2, nz);

[x, y, z] = ndgrid(vx, vy, vz); % ndgrid of position vectors

ux = 5;

u = ux  .* ones(size(x));
v = zeros(size(x));
w = zeros(size(x));


%%

% quiver3(x, y, z, u, v, w)


%% Have the co-ordinate system in polar form

r = sqrt(x.^2 + y.^2 + z.^2);
theta = atan(z./x);
phi = atan(y./x);

% [r, theta, phi] = ndgrid(r, theta, phi);

vr = sin(theta) .* cos(phi) .* u + sin(theta) .* sin(phi) .* v + cos(theta) .* w;
vtheta = cos(theta) .* cos(phi) .* u + cos(theta) .* sin(phi) .* v - sin(theta) .* w;
vphi = -sin(phi) .* u + cos(phi) .* v;

range = sqrt(vx(end).^2 + vy(end).^2);

% z_ = 0;
% 
% [~, indz] = find(abs(vz - z_) < 1e2);
% x_ = x(:, :, indz(end));
% y_ = y(:, :, indz(end));

idx = (abs(theta(:) - pi/2) < 1e1) & (abs(r(:) - 5e2) < 1e-3);

idx = reshape(idx, size(theta));

% [idxx idxy idxz] = ind2sub(size(theta),c);

figure; plot(x(idx == 1), y(idx == 1), 'o');


figure; plot(vr(idx == 1), 'o');