clear;
close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1.8;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
% Omega_rpm = linspace(1, 60, 4);
Omega_rpm = 60;
% Omega_rpm = 1;

lambda = 0.03;
%% 

for i = 1:length(Omega_rpm)
   
    
BW = BW_deg * pi/180;
Omega(i) = Omega_rpm(i) * 2*pi/60; % rotation speed in rad/s
Td = BW/Omega(i);
hs = round(Td/PRT);

hs_(i) = hs;

% if mod(hs, 2) == 0
%     hs = hs + 1;
% end

phi_0 = phi_0_deg * pi/180;

sec = round((n_rot*2*pi)/BW);

N = sec * hs;
t = 0:PRT:(N - 1)*PRT;
% phi_axis = phi_0:BW:phi_0+Omega*t(end);

phi_axis_1 = phi_0:BW:phi_0+n_rot.*2*pi;
phi_axis = mean([phi_axis_1(1:end-1); phi_axis_1(2:end)]);

% phi_axis = zeros(size(phi_axis));

%% variables for signal model

beta_wind_deg = 0;

beta_wind = beta_wind_deg .* pi/180;

v_amb = lambda/(4 * PRT);

% mu_mean = linspace(-v_amb, v_amb, 100);
mu_mean = 5;
for l = 1:length(mu_mean)

mu = normrnd(mu_mean(l), 1, [1 N]);
% mu = normrnd(4, 1, [1 N]);
n_sig = sqrt(0.01);

% mu = 3;


[s, SNR] = TD_generator(mu, lambda, beta_wind, phi_0, Omega(i), t, n_sig);

N = length(t);



if mod(hs, 2) == 0
    vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
%      vel_axis_hs = linspace(-v_amb, v_amb, hs+1);
else
    vel_axis_hs = linspace(-v_amb, v_amb, hs);
end

[S1, F1, Ti1, P1] = spectrogram(s, hs, 0, hs , N*1/60);


% S1 = S1./max(S1(:));
S1_norm = 1./sqrt(hs) .* fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));


txt = ['SNR = ', num2str(db(SNR)/2), ' dB , \Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];



figure;
imagesc(vel_axis_hs, phi_axis*180/pi, S1_norm_db); shading flat;
colormap('jet');
colorbar;


xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title(txt, 'FontSize', 12, 'FontWeight', 'bold');




end
end