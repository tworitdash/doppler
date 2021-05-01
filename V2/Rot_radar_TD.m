clear;
close all;
%% RADAR constants and wind direction
SNR_db = 30;
SNR = 10^(SNR_db/10);

BW_deg = 1.8;
BW = BW_deg * pi/180;


beta_wind = eps; % wind direction
mu = 5;
sigma = 0.2;



PRT = 1e-3;
lambda = 3e-2;
n = 2^10;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity


Omega_rpm = 1;


Phi_0_deg = 0;
Phi_end_deg = 360;
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle
phi_end = Phi_end_deg * pi/180; % End of azimuth angle

Phi = phi_0:BW:phi_end;

    

Omega = Omega_rpm .* 2 * pi ./ 60;

T = BW/Omega;

time_axis = eps:PRT:T;
            

hits_scan_ = length(time_axis); % length of time axis

% hits_scan_ = 2^(nextpow2(hits_scan) - 1);
          
delta_v = lambda/(2*hits_scan_*PRT);
            
            

for k = 1:length(Phi)-1
      beta_scan = beta_wind - linspace(Phi(k), Phi(k + 1), hits_scan_);
      beta_scan_hd = beta_wind - linspace(Phi(k), Phi(k + 1), n);
      [data, data_f] = DS_simulatorV2(SNR, 1, mu, sigma, n, v_amb, hits_scan_);
      
      sig(k, :) = (abs(squeeze(data))...
                 .* exp(1j .* unwrap(angle(squeeze(data))) .* cos(beta_scan)));
      
end
            
   
signal = sig(1, :);
I = real(signal); 
Q = imag(signal);
dt = PRT;

figure; plot(time_axis, I, '-*'); hold on; plot(time_axis, Q, '-o');

for i = 1:size(sig, 2) - 1
   Num_int(i) = (Q(i + 1) .* I(i) - I(i + 1) .* Q(i));
   Denum_int(i) = (I(i).^2 + Q(i).^2);
end

omega_mean = 1./dt .* sum(Num_int)./sum(Denum_int);
v_mean = omega_mean .* lambda ./ 2;

