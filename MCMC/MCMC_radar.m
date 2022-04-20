%% Radar Data Generator
clear;
close all;

markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('../Random_Fields/colors.mat');
colors = colors.color;

%% Define the 3D resolution cell in terms of r, theta and phi with scatterers positions
% Cell parameters


r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 1.8*pi/180;
theta_0 = 0*pi/180;
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

% figure; scatter3(x0, y0, z0, 50, r0); colormap('copper'); colorbar; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 

%% Signal model

D_min = 0.6;
D_max = 4;
D0 = 1;
N0 = 8e3;
u_mean = 1;
v_mean = 0;
u_sigma = 1;
v_sigma = 0;
SNR_db = 30;
n_rot = 15;
NSweep = 10;
NSec = 100;
Nt = n_rot * NSweep * NSec;
% Nt = 10;
dt = 1e-3;
lambda = 3e-2;

[Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma, sigma_n] = ...
    Signal_modelvmac(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);

% figure; plot(real(Z_model)); hold on; plot(imag(Z_model));
figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));


Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] = Spec(Z_model_re(1:NSweep), NSweep, dt, lambda, SNR_db, 1, 1, 7);

%% Signal properties

Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

% figure; plot(real(Z_avail)); hold on; plot(imag(Z_avail));
figure; histogram(real(Z_avail)); figure; histogram(imag(Z_avail));

% Z_avail_vec is the signal available (measurements)


%% The model for MCMC

u_mean_obs = mu_obs; 
u_sigma_obs = sigma_obs;

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep];
end

t_avail = reshape(t, [length(Z_avail_vec) 1]);

[accepted_umean, rejected_umean, accepted_usigma, rejected_usigma] = MHradar([u_mean_obs, u_sigma_obs], 200000, (Z_avail_vec), t_avail, x0, y0, N, sigma_n);

figure; plot(rejected_usigma); hold on; plot(accepted_usigma);

hs = mean(accepted_usigma);

figure; plot(rejected_umean); hold on; plot(accepted_umean);

hm = mean(accepted_umean);

figure; histogram(accepted_umean); figure; histogram(accepted_usigma);

%% Reconstruction 

U_re = normrnd(hm, hs, [1 N]);

t_re = eps:dt:(Nt - 1)*dt;

for l = 2:length(t_re)
    r_re = sqrt((x0 + U_re .* t_re(l)).^2 + y0.^2);
    Z_re(l) = sum(exp(1j .* 4 .* pi/lambda .* r_re)) + sigma_n.^2 .* randn;
end

% Spec(Z_re, Nt, dt, lambda, SNR_db, 1, 1, 1)

figure(1); hold on; histogram(real(Z_re)); figure(2); hold on; histogram(imag(Z_re));
