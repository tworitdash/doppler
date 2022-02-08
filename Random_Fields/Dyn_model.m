%% Signal model with HD

clear;
close all;

%% 

markers = load('../mono/markers.mat');
markers = markers.markers;

N = 100;
r0 = 100;
dr = 50;
ph0 = 0;
n_bw = 2; % number of beamwidths needed
dph = n_bw*1.8*pi/180;
n_rot = 32;
N_Sweep = 5;
N_sec = 200;
Nt = n_rot*N_sec*N_Sweep;
SNR_db = 30;
lambda = 0.03;
dt = 1e-3;
u_mean = 0;
u_sigma = 1;
v_mean = 2;
v_sigma = 0;

%% Generating measurements with only Nt samples per rotation

[z_model] = Zmodel(N, r0, dr, ph0, dph, Nt, SNR_db, lambda, dt, u_mean, u_sigma, v_mean, v_sigma); 

[ZFFT, PT, mu, sigma, vel_axis, dv] = Spec(z_model, Nt, dt, lambda, SNR_db, 1, 1, 5);

[zi, z] = Zavail(z_model, n_rot, n_bw*N_Sweep, N_sec/n_bw);

% [~,~, mu_avail, sigma_avail, vel_axis_avail, dv_avail] = Spec(zi(2, :), N_Sweep, dt, lambda, SNR_db, 1, 1, 1);

%% Processing with model and measurements
%% Initial data processing to find rough values of Doppler moments

n_bw_p = 2; % number of beamwidth you want to process
N_Sweep_tot = n_bw_p * N_Sweep;
N_Sweep_p = 5; % Number of sweeps you want to integrate
N_Sweep_Skip = round(N_Sweep_tot/N_Sweep_p); % Number of Sweeps that will be skipped

n_rot_p = 8; % Number of rotations you want to average
dtp = N_Sweep_Skip * dt; % New pulse repition time

for i = 1:n_rot_p

data0 = zi(i, 1:N_Sweep_Skip:N_Sweep_tot);

if i == 1
    [~,~, mu_avail(i), sigma_avail(i), vel_axis_avail, dv_avail] = Spec(data0, N_Sweep_p, dtp, lambda, SNR_db, 1, 1, 6);
else
    [~,~, mu_avail(i), sigma_avail(i), vel_axis_avail, dv_avail] = Spec(data0, N_Sweep_p, dtp, lambda, SNR_db, 0, 1, 2);
end

end

mu_emp = mean(mu_avail);
sigma_emp = mean(sigma_avail);

%% Generating data for the time gaps between rotations

% Spatial resolution wanted

Np = 100;
r0p = 100;
drp = 50;
ph0p = 0;
n_bwp = 1; % number of beamwidths needed
dphp = n_bw*1.8*pi/180;

N_Sweep_avail = n_bwp*N_Sweep;
N_sec_avail = N_sec/n_bwp;

N_Sweepp = 5; % Number of sweeps you need 

N_Sweep_Skipp = round(N_Sweep_avail/N_Sweepp); % Number of Sweeps that will be skipped

N_gap = (N_sec_avail - 1)*N_Sweep_avail;
dtpp = 1e-3;

mu_model = mu_emp;
sigma_model = sigma_emp; Spec(data0, N_Sweep_p, dtp, lambda, SNR_db, 1, 1, 6);

[z_modelp] = Zmodel(N, r0p, drp, ph0p, dphp, N_gap, SNR_db, lambda, dtpp, mu_model, sigma_model, 0, 0); 

%% Data with rotation number n_rot_p and n_rot_p+1 with simulated data in between them

n_rotpp = 2; % Mumber of rotations you want to integrate
z_nrotp = zi(n_rot_p, 1:N_Sweep_Skipp:N_Sweepp);

z_hybrid = z_nrotp;
% 
for l = n_rot_p+1:n_rot_p+n_rotpp

z_nrotp1 = zi(l, 1:N_Sweep_Skipp:N_Sweepp);
z_hybrid = [z_hybrid z_hybrid(end).*z_modelp z_nrotp1];
% z_hybrid = [z_hybrid  z_nrotp1];
% 
% z_nrotp = zi(n_rot_p, 1:N_Sweep_Skipp:N_Sweepp);
% z_nrotp1 = zi(n_rot_p+1, 1:N_Sweep_Skipp:N_Sweepp);
% z_nrotp2 = zi(n_rot_p+2, 1:N_Sweep_Skipp:N_Sweepp);

end

% z_hybrid = hann(1, length(z_hybrid)) .* z_hybrid;

% z_hybrid = [z_nrotp z_nrotp(end).*z_modelp z_nrotp1 z_nrotp1(end).*z_modelp z_nrotp2];

[ZFFT, PT, mu_re, sigma_re, vel_axis_re, dv_re] = Spec(z_hybrid, length(z_hybrid), dtpp, lambda, SNR_db, 1, 1, 3);
