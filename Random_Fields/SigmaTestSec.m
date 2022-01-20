%% Signal model with HD
clear;
% close all;

markers = load('../mono/markers.mat');
markers = markers.markers;
U_sigma = linspace(eps, 7.5, 100);
MC = 32;
for k = 1:MC
for i = 1:length(U_sigma)

N = 100;
r0 = 100;
dr = 50;
ph0 = 0;
dph = 2*1.8*pi/180;
n_rot = 1;
N_Sweep = 5;
N_sec = 200;
Nt = n_rot*N_sec*N_Sweep;
SNR_db = 30;
lambda = 0.03;
dt = 2*1e-3;
u_mean = 0;
u_sigma = U_sigma(i);
v_mean = 0;
v_sigma = 0;

[z_model] = Zmodel(N, r0, dr, ph0, dph, Nt, SNR_db, lambda, dt, u_mean, u_sigma, v_mean, v_sigma); 

[ZFFT, PT, mu(i, k), sigma(i, k), vel_axis, dv, v_amb] = Spec(z_model, Nt, dt, lambda, SNR_db, 0, 2);

%% Decide the samples available and club them

[zi, z] = Zavail(z_model, n_rot, N_Sweep, N_sec);

[~,~, mu_avail(i, k), sigma_avail(i, k), vel_axis_avail, dv_avail] = Spec(zi(1, :), N_Sweep, dt, lambda, SNR_db, 0, 1);

end
end

Mu_True_Avg = mean(mu, 2);
Sigma_True_Avg = mean(sigma, 2);
Mu_Avail_Avg = mean(mu_avail, 2);
Sigma_Avail_Avg = mean(sigma_avail, 2);

%% Plot with small number of samples

txt = ['Mean Doppler Velocity'];
dtext = ['Nt = ', num2str(N_Sweep), ', dv = ', num2str(dv_avail), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]', ', d\phi = ', num2str(dph*180/pi), '^{\circ}' ];
xl = 'True Spectrum Width \sigma_{True}';

f = figure(1001); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['Mean Doppler Velocity [m/s]'];

marker = markers(1);
plott2(U_sigma, Mu_Avail_Avg, xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);

%% Plot with large number of samples


dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]', ', d\phi = ', num2str(dph*180/pi), '^{\circ}' ];
marker = markers(2);
plott2(U_sigma, Mu_True_Avg, xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);

%% Plot with small number of samples


txt = ['Doppler spectrum width'];
dtext = ['Nt = ', num2str(N_Sweep), ', dv = ', num2str(dv_avail), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]', ', d\phi = ', num2str(dph*180/pi), '^{\circ}' ];
% xl = 'True Normalized Spectrum Width \sigma_{True}/v_{amb}';
xl = 'True Spectrum Width \sigma_{True}';
f = figure(1000); hold on; f.Position = [10 10 1000 1000];
color = 'k';
% yl =  ['Doppler Spectrum Width Normalized \sigma_{retrieved}/v_amb [m/s]'];
yl =  ['Doppler Spectrum Width [m/s]'];
marker = markers(1);
plott2(U_sigma, Sigma_Avail_Avg, xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);

%% Plot with large number of samples

dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]', ', d\phi = ', num2str(dph*180/pi), '^{\circ}' ];
marker = markers(2);

plott2(U_sigma, Sigma_True_Avg, xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);

%% 