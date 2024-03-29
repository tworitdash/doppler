clear;
% close all;

markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

%% Define the 3D resolution cell in terms of r, theta and phi with scatterers positions
% Cell parameters


r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 1*pi/180;
theta_0 = 0*pi/180;
dth = 2*pi/180;

% Scatterers positions

N = 2000;
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
% u_mean_ = linspace(-7.5, 7.5, 100);
u_sigma_ = linspace(0.5, 7.5, 100);
% N_Sweep = [300 15 9 5];
N_Sweep = 5;
Omega = [1 20 40 60];

for k = 1:length(N_Sweep)

for i = 1:length(u_sigma_)
u_mean = 0;
v_mean = 0;
u_sigma = u_sigma_(i);
v_sigma = 0;
SNR_db = 30;
n_rot = 32;
NSweep = N_Sweep(k);
NSec = 1;
Nt = n_rot * NSweep * NSec;
% Nt = 10;
dt = 1e-3;
lambda = 3e-2;

[Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma] = ...
    Signal_modelvmac(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);


%% processing

Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);

for l = 1:n_rot

[ZFFT_, PT_(l), mu__(l), sigma__(l), vel_axis_, dv_] = Spec(Z_model_re(1:NSweep, l).', NSweep, dt, lambda, SNR_db, 0, 1, 7);

end

% txt = ['Doppler Spectrum'];
% % dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
% dtext = ['N_{sweeps} = ', num2str(NSweep), ', dv = ', num2str(dv_), ' [m/s]'];
% xl = 'velocity [m/s]';
% 
% % f = figure(108); hold on; f.Position = [10 10 1000 1000];
% 
% %zl =  ['Angle [{\circ}]'];
% yl =  ['Spectrum [dB]'];
% color = colors(2).c;
% marker = markers(k);
% plott2(vel_axis_, db(abs(ZFFT_)), xl, yl, txt, 2, dtext, color, marker); 

mu_(k, i) = mean(mu__);
sigma_(k, i) = mean(sigma__);
end

txt = ['SNR = ', num2str(SNR_db), ' [dB]'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ];
dtext = ['BW = ', num2str(dph*180/pi), '^{\circ}'];
xl = '\sigma True [m.sec^{-1}]';

f = figure(108); hold on; f.Position = [10 10 1000 1000];
color = colors(4).c;
yl =  ['\mu retrieved [m.sec^{-1}]'];

marker = markers(k);
plott2(u_sigma_, mu_(k, :), xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);

txt = ['SNR = ', num2str(SNR_db), ' [dB]'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ];
dtext = ['BW = ', num2str(dph*180/pi), '^{\circ}'];
xl = '\sigma True [m.sec^{-1}]';

f = figure(109); hold on; f.Position = [10 10 1000 1000];
color = colors(4).c;
yl =  ['\sigma retrieved [m.sec^{-1}]'];

marker = markers(k);
plott2(u_sigma_, sigma_(k, :), xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);

end