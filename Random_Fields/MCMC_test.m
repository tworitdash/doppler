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

N = 1000;
r = r0-dr/2 + dr* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);
th = theta_0 - dth/2 + dth .* rand(1, N);

% positions of these scatteres in Cartesian domain

x0 = r .* cos(th) .* cos(ph);
y0 = r .* cos(th) .* sin(ph);
z0 = r .* sin(th); 

r0 = sqrt(x0.^2 + y0.^2 + z0.^2);
figure; scatter3(x0, y0, z0, 50, r0); colormap('copper'); colorbar; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 

%% Signal model

D_min = 0.6;
D_max = 4;
D0 = 1;
N0 = 8e3;
N = 1000; 
u_mean = 1;
v_mean = 1;
u_sigma = 2;
v_sigma = 2;
SNR_db = 30;
n_rot = 10;
NSweep = 5;
NSec = 200;
Nt = n_rot * NSweep * NSec;
dt = 1e-3;
lambda = 3e-2;

[Z_model, v_amb, Vtmean, Vtspread, vr, u, v, Wt, D] = ...
    Signal_model(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);

[ZFFT, PT, mu, sigma, vel_axis, dv] = Spec(Z_model(1:128), 128, dt, lambda, SNR_db, 1, 1, 5);


[p_u, u_x] = hist(u);
[p_v, v_x] = hist(u);
[p_W, W_x] = hist(Wt);
P_uvW = p_u .* p_v .* p_W;

figure; scatter3(u_x, v_x, W_x, 50, P_uvW); colormap('copper'); colorbar; xlabel('u [m/s]'); ylabel('v [m/s]'); zlabel('Wt [m/s]'); 


%% Signal available 
Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);

Z_model_re_re = Z_model_re(1:NSweep, :);

Z_model_meas = reshape(Z_model_re_re, [1 n_rot * NSweep]);

%% Available Spectrum
n_rot_a = 1;
Z_model_meas_a = Z_model_meas(1:n_rot_a * NSweep);

[ZFFTm, PTm, mum, sigmam, vel_axism, dvm] = Spec(Z_model_meas_a, length(Z_model_meas_a), dt, lambda, SNR_db, 1, 1, 6);

%% MCMC