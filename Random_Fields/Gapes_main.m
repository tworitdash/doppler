clear;
close all;

%% 

markers = load('../mono/markers.mat');
markers = markers.markers;

N = 10000;
r0 = 100;
dr = 50;
ph0 = 0;
n_bw = 1; % number of beamwidths needed
dph = n_bw*1.8*pi/180;
n_rot = 2;
N_Sweep = 10;
N_sec = 100;
Nt = n_rot*N_sec*N_Sweep;
SNR_db = 30;
lambda = 0.03;
dt = 1e-3;
u_mean = 0;
u_sigma = 2;
v_mean = 0;
v_sigma = 0;

%% Generating measurements with only Nt samples per rotation

[z_model] = Zmodel(N, r0, dr, ph0, dph, Nt, SNR_db, lambda, dt, u_mean, u_sigma, v_mean, v_sigma); 

[ZFFT, PT, mu, sigma, vel_axis, dv] = Spec(z_model, Nt, dt, lambda, SNR_db, 1, 1, 5);

[zi, z, z_with_gap] = Zavail(z_model, n_rot, n_bw*N_Sweep, N_sec/n_bw);
%% Test Gapes

sig = [z_with_gap(1, :) z_with_gap(2, :)];
M = 4;
N_omg = length(z_model);

for i = 1:n_rot
    Nvec((i - 1)*N_Sweep+1:i*N_Sweep) = N_sec*N_Sweep*(i - 1)+1:N_sec*N_Sweep*(i - 1)+N_Sweep;
end

[H, beta, sig_out] = GAPES(sig, Nvec, M, N_omg);

[ZFFTa, PTa, mua, sigmaa, vel_axisa, dva] = Spec(sig(1:N_Sweep), N_Sweep, dt, lambda, SNR_db, 1, 1, 7);

[ZFFTre, PTre, mure, sigmare, vel_axisre, dvre]  = Spec(sig_out.', length(sig_out), dt, lambda, SNR_db, 1, 1, 8);
