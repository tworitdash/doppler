%% Weather simulator

clear; 
close all;

input_info.r0 = 100;
input_info.phi0 = 0*pi/180;
input_info.dr = 0;
input_info.dph = 0; % 0.001*pi/180;

input_info.spatial_dist.type = 1;
input_info.plot_geo = 1;
input_info.NScatters = 2;

input_info.spatial_dist.lamr = 10;
input_info.spatial_dist.lamp = 20;

input_info.RADAR.dT = 1e-3;

input_info.RADAR.dTjittervec = 0; %input_info.RADAR.dT/3; % linspace(0, 4*input_info.RADAR.dT, 5);
input_info.N_pulse = 5;

input_info.Ngap_avg = 95;

input_info.sig_gap = 0; %2 * input_info.N_pulse .* linspace(0, 2, 3);

for k = 1:length( input_info.sig_gap )
    
input_info.RADAR.dTjitter = 0; % input_info.RADAR.dTjittervec(k); % input_info.RADAR.dT/2;
input_info.RADAR.lambda = 3e-2;


input_info.N_rot         = 32;

input_info.N_gap         = [0 randi([input_info.Ngap_avg-input_info.sig_gap(k)/2 input_info.Ngap_avg+input_info.sig_gap(k)/2], 1, input_info.N_rot)];
input_info.velocity.u.mu = 3;
input_info.velocity.v.mu = 0;
input_info.velocity.u.sigma = 1;
input_info.velocity.v.sigma = 0;

input_info.SNR = 20;
input_info.velocity.type = 2;
input_info.Doppler_plot = 0;
% input_info.Ngt = 128;
input_info.Ngt = [];


[out] = SimRad(input_info); % Signal generator

[ZFFT, PT, mu, sigma, vel_axis, dv, v_amb] = Spec(out.Z_avail_vec, length(out.Z_avail_vec), out.dTavg, input_info.RADAR.lambda, input_info.SNR, 1, 1, 1);

end
%% 

