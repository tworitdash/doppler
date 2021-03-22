close all;
clear;
%% Created grid of space
[x, y, z] = ndgrid(linspace(-2, 2, 100), 0, linspace(0, 2, 100));

%% Point out where vortices are


x0 = 1; % Vortex x co-ordinates
y0 = 0; % Vortex y co-ordinates
z0 = 1; % Vortex z co-ordinates

rx0 = 1; % Vortex x co-ordinates
ry0 = 0; % Vortex y co-ordinates
rz0 = 1; % Vortex z co-ordinates

%% scaled co-ordinates

x_ = (x - x0)/rx0;
% y_ = (y - y0)/ry0; 
y_ = y;
z_ = (z - z0)/rz0;


theta_plus = 0;
phi_plus = 0;

r_ = sqrt(x_.^2 + y_.^2 + z_.^2);
theta_ = pi/2 - atan(z_./r_) + theta_plus;
phi_ = atan(y_./x_) + phi_plus;


A0 = 1;
sigma = 1;

A_ = A0  .* exp(-(r_ - 1)/(2.* sigma.^2));

u = A_ .* cos(theta_) .* cos(phi_);
v = A_ .* cos(theta_) .* sin(phi_);
w = -A_ .* sin(theta_);


e = 1:5:size(x, 1);
figure; quiver3(x(e, 1, e), y(e, 1, e), z(e, 1, e), u(e, 1, e), v(e, 1, e), w(e, 1, e), 'AutoScale','off');
xlabel('x');
ylabel('y');
zlabel('z');