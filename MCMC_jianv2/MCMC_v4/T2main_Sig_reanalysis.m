clear; 
% close all;

input_info.r0 = 100;
input_info.phi0 = 0*pi/180;
input_info.dr = 50;
input_info.dph = 2*pi/180; % 0.001*pi/180;

input_info.spatial_dist.type = 3;
input_info.plot_geo = 0;
input_info.NScatters = 10000;

input_info.spatial_dist.lamr = 1;
input_info.spatial_dist.lamp = 2;

input_info.RADAR.dT = 1e-3;

input_info.RADAR.dTjittervec = 0; %input_info.RADAR.dT/3; % linspace(0, 4*input_info.RADAR.dT, 5);


input_info.Ngap_avg = 0;

input_info.sig_gap = 0; %32 * input_info.N_pulse; %2 * input_info.N_pulse .* linspace(0, 2, 3);1
MC = 128;
%%

for k = 1:length( input_info.sig_gap )
    
input_info.RADAR.dTjitter = 0; % input_info.RADAR.dTjittervec(k); % input_info.RADAR.dT/2;
input_info.RADAR.lambda   = 3e-2;


input_info.N_rot          = 1;

input_info.N_gap          = [0 randi([input_info.Ngap_avg-input_info.sig_gap(k)/2 input_info.Ngap_avg+input_info.sig_gap(k)/2], 1, input_info.N_rot)];
input_info.velocity.u.mu  = 7.5;
input_info.velocity.v.mu  = 0;

input_info.v_amb = input_info.RADAR.lambda/(2.*input_info.RADAR.dT);
input_info.velocity.v.sigma = 0;


input_info.velocity.u.sigmavec = linspace(0, input_info.v_amb/2, 1000);
input_info.N_pulsevec = 5:1:128;



input_info.SNR = 30;
input_info.velocity.type = 2;
input_info.Doppler_plot = 0;
input_info.Ngt = 128;
% input_info.Ngt = [];
input_info.Doppler_plot = 1;
input_info.vel_amb = 0;
% input_info.plot_geo = 1;

% rng(1, 'twister');
for s = 1:length(input_info.velocity.u.sigmavec)
    
for n = 1:length(input_info.N_pulsevec)
    
input_info.velocity.u.sigma = input_info.velocity.u.sigmavec(s);
input_info.N_pulse = input_info.N_pulsevec(n);


parfor mc = 1:MC
    [out] = SimRad(input_info); % Signal generator
    mu_re(s, n, mc) = out.mu_obszp;
    sig_re(s, n, mc) = out.sigma_obszp;
end

Mu_error(s, n) = sqrt( 1/MC .* sum( ( squeeze(mu_re(s, n, :)) - input_info.velocity.u.mu ).^2 ));
Sig_error(s, n) = sqrt( 1/MC .* sum( ( squeeze(sig_re(s, n, :)) - input_info.velocity.u.sigma ).^2 ));

Mu(s, n) = mean( squeeze(mu_re(s, n, :)) );
Sigma(s, n) = mean( squeeze(sig_re(s, n, :)) );

% dv_re(s, n) = out.dv_obs;
% figure(1001); hold on; plot(real(out.z)); hold on; plot(imag(out.z));

end
end

end
%%
x = input_info.N_pulsevec;
y = input_info.velocity.u.sigmavec;
z = Mu_error;

xl = [' Number of Pulses ' ];
yl = [' True velocity width u_{sigma}' ];
zl = ['Error in mean velocity estimation u_{\mu_{re}}' ];
tl = ['Error in mean velocity estimation u_{\mu_{re}}' ];
surplot_pcolor(x, y, z, xl, yl, zl, tl)

x = input_info.N_pulsevec;
y = input_info.velocity.u.sigmavec;
z = Sig_error;

xl = [' Number of Pulses ' ];
yl = [' True velocity width u_{sigma}' ];
zl = ['Error in velocity width estimation u_{\sigma_{re}}' ];
tl = ['Error in velocity width estimation u_{\sigma_{re}}' ];
surplot_pcolor(x, y, z, xl, yl, zl, tl)

