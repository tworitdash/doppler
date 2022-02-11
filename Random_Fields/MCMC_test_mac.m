clear;
close all;

markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

%% Define the 3D resolution cell in terms of r, theta and phi with scatterers positions
% Cell parameters

r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 1.8*pi/180;
theta_0 = 30*pi/180;
dth = 2*pi/180;

% Scatterers positions

N = 200;
r = r0-dr/2 + dr* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);
th = theta_0 - dth/2 + dth .* rand(1, N);

% positions of these scatteres in Cartesian domain

x0 = r .* cos(th) .* cos(ph);
y0 = r .* cos(th) .* sin(ph);
z0 = r .* sin(th); 

% r0 = sqrt(x0.^2 + y0.^2 + z0.^2);
r0 = sqrt(x0.^2 + y0.^2); % + z0.^2);

figure; scatter3(x0, y0, z0, 50, r0); colormap('copper'); colorbar; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 

%% Signal model

D_min = 0.6;
D_max = 4;
D0 = 1;
N0 = 8e3;
u_mean = 1;
v_mean = 1;
u_sigma = 2;
v_sigma = 2;
SNR_db = 30;
n_rot = 5;
NSweep = 10;
NSec = 100;
Nt = n_rot * NSweep * NSec;
% Nt = 10;
dt = 1e-3;
lambda = 3e-2;

[Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma] = ...
    Signal_modelvmac(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);

figure; plot(real(Z_model)); hold on; plot(imag(Z_model));
figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));

%% Signal available 
Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);

Spec(Z_model_re(1:NSweep), NSweep, dt, lambda, SNR_db, 1, 1, 7);

%% Signal properties

Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

figure; plot(real(Z_avail)); hold on; plot(imag(Z_avail));
figure; histogram(real(Z_avail)); figure; histogram(imag(Z_avail));


%% Prior


u_mean_prior = 0.1;
v_mean_prior = 0;
u_sigma_prior = 1;
v_sigma_prior = 1;

W_mean_prior = Vtmean;
W_sigma_prior = Vtspread;

likelihood_current = cumprod(normpdf(real(Z_avail_vec), u_mean_prior, u_sigma_prior));

