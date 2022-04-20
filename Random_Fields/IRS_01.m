clear;
% close all;

markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

%% Define the 3D resolution cell in terms of r, theta and phi with scatterers positions
% Cell parameters
ph_ = linspace(eps, 360, 360);
% ph_ = 0;

for i = 1:length(ph_)

r0 = 100;
dr = 50;
ph0 = ph_(i)*pi/180;
dph = 1*pi/180;
theta_0 = 0*pi/180;
dth = 2*pi/180;

% Scatterers positions

N = 2;
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
u_mean = 5;
v_mean = 0;
u_sigma = 1;
v_sigma = 0;
SNR_db = 30;
n_rot = 1;
NSweep = 5;
NSec = 1;
Nt = n_rot * NSweep * NSec;
% Nt = 10;
dt = 1e-3;
lambda = 3e-2;

[Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma] = ...
    Signal_modelvmac(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);


%% processing

Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);

[ZFFT_(i, :), PT_(i), mu_(i), sigma_(i), vel_axis_, dv_] = Spec(Z_model_re(1:NSweep).', NSweep, dt, lambda, SNR_db, 0, 1, 7);

end

tl = ['Doppler Spectrum'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
dtext = ['N_Sweep = ', num2str(NSweep), ', dv = ', num2str(dv_), ' [m/s]', ', \sigma_{re} = ', num2str(sigma_), ' [m/s]'];
xl = 'velocity [m/s]';

% f = figure(108); hold on; f.Position = [10 10 1000 1000];

yl =  ['Angle [{\circ}]'];
zl =  ['Spectrum [dB]'];

surplot_surf(vel_axis_, ph_, db(abs(ZFFT_)), xl, yl, zl, tl);


