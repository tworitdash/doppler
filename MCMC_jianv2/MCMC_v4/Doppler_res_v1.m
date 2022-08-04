clear;
markers = load('markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

%% Weather simulator

 
% close all;

input_info.r0 = 100;
input_info.phi0 = 0*pi/180;
input_info.dr = 0;
input_info.dph = 0*pi/180; % 0.001*pi/180;

input_info.spatial_dist.type = 1;
input_info.plot_geo = 0;
input_info.NScatters = 200000;

input_info.spatial_dist.lamr = 10;
input_info.spatial_dist.lamp = 20;

input_info.RADAR.dT = 1e-3;

input_info.RADAR.dTjittervec = 0; %input_info.RADAR.dT/3; % linspace(0, 4*input_info.RADAR.dT, 5);

input_info.N_pulsevec = 1024; % 5:1:128;

for n = 1:length(input_info.N_pulsevec)

input_info.N_pulse = input_info.N_pulsevec(n);

input_info.Ngap_avg = 0;

input_info.sig_gap = 0; %32 * input_info.N_pulse; %2 * input_info.N_pulse .* linspace(0, 2, 3);1

%%


disp(n);

for k = 1:length( input_info.sig_gap )
    
input_info.RADAR.dTjitter = 0; % input_info.RADAR.dTjittervec(k); % input_info.RADAR.dT/2;
input_info.RADAR.lambda   = 3e-2;


input_info.N_rot          = 1;

input_info.N_gap          = [0 randi([input_info.Ngap_avg-input_info.sig_gap(k)/2 input_info.Ngap_avg+input_info.sig_gap(k)/2], 1, input_info.N_rot)];
input_info.velocity.u.mu  = 7.5;
input_info.velocity.v.mu  = 0;


input_info.velocity.u.sigmavec = linspace(eps, 4, 100);

for i = 1:length(input_info.velocity.u.sigmavec)

disp(i);
    
input_info.velocity.u.sigma = input_info.velocity.u.sigmavec(i);
input_info.velocity.v.sigma = 0;

input_info.SNR = 3000;
input_info.velocity.type = 2;
input_info.Doppler_plot = 0;
input_info.Ngt = 128;
% input_info.Ngt = [];
input_info.Doppler_plot = 1;
input_info.vel_amb = 0;
% input_info.plot_geo = 1;

% rng(1, 'twister');


    
[out] = SimRadGauss(input_info); % Signal generator

Mu(n, i) = out.mu_obs;
Sigma(n, i) = out.sigma_obs;

end

end
end
%% 

txt = ['\sigma_{True} vs \sigma_{re}'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]'];

dtext = ['Nt = ', num2str(input_info.N_pulse)];
% ', dv = ', num2str(dv_obs), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]', ...
%    ', \mu_{True} = ', num2str(mu_true), ' [m/s]', ', \sigma_{True} = ', num2str(sigma_true), ' [m/s]'];

xl = '\sigma_{True} [m/s]';

f = figure(1000); hold on; f.Position = [10 10 1000 1000];
color = colors(3).c;
yl =  ['\sigma_{re}'];

marker = markers(1);
plott2(input_info.velocity.u.sigmavec, Sigma, xl, yl, txt, 2, dtext, color, marker);

marker = markers(2);

plott2(out.dv_obs .* ones(1, length(input_info.velocity.u.sigmavec)), ...
    Sigma, xl, yl, txt, 2, dtext, color, marker);

% plott2(x, y, xl, yl, tl, lw, dtext, c, m)
%% Error Plots
% 
% x = input_info.N_pulsevec;
% y = input_info.velocity.u.sigmavec;
% z = Mu_error;
% 
% xl = [' Number of Pulses ' ];
% yl = [' True velocity width u_{sigma}' ];
% zl = ['Error in mean velocity estimation u_{\mu_{re}}' ];
% tl = ['Error in mean velocity estimation u_{\mu_{re}}' ];
% surplot_pcolor(x, y, z.', xl, yl, zl, tl)
% 
% 
% x = input_info.N_pulsevec;
% y = input_info.velocity.u.sigmavec;
% z = Sigma_error;
% 
% xl = [' Number of Pulses ' ];
% yl = [' True velocity width u_{sigma}' ];
% zl = ['Error in velocity width estimation u_{\sigma_{re}}' ];
% tl = ['Error in velocity width estimation u_{\sigma_{re}}' ];
% surplot_pcolor(x, y, z.', xl, yl, zl, tl)