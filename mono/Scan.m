%% Radar simulator on an observation rim [Range is constant]
clear; 
close all;
%% Radar paramaters
lambda = 0.03;
PRT = 1e-3; % Pulse repitition time
v_amb = lambda/(4 * PRT);

Omega_rpm = 1; % Scanning rate of radar in azimuth in rpm
Omega = Omega_rpm * 2 * pi / 60;

BW_deg = 1.8; % BeamWidth in degree
BW = BW_deg * pi/180;

sec = round(2*pi/BW); % number of sectors

Td = BW/Omega; 
hs = round(Td/PRT);
N = hs * sec;
phi_0 = 0; % start angle of the scan

%% Wind Parameters

u = 5; % wind speed mean in space
sigma = 1; % wind speed standard deviation in space

U = normrnd(u, sigma, [1 N]); % Gaussian distrbution of wind speed over the circular rim
phi_wind = 0;

phi_axis_1 = phi_0:BW:phi_0+2*pi;
phi_axis = mean([phi_axis_1(1:end-1); phi_axis_1(2:end)]);

%% Time axis/ Noise/ Radar simulator 

t = 0:PRT:(N - 1)*PRT;
n_sig = sqrt(1e-3);
ph = 4*pi/lambda .* U .* cos(phi_wind - Omega .* t - phi_0) .* t;
s = exp(1j * ph) + n_sig .* randn(size(t)); % radar signal
% [s, SNR] = TD_generator(U, lambda, phi_wind, phi_0, Omega, t, n_sig);
% SNR = sum(abs(s).^2)/(length(t) .* n_sig^2);

%% Retrieval 

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


txt = ['SNR = ', num2str(db(SNR)/2), ' dB , \Omega = ', num2str(Omega_rpm), ' [rpm]'];



figure;
imagesc(vel_axis_hs, phi_axis*180/pi, S1_norm_db); shading flat;
colormap('jet');
colorbar;


xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title(txt, 'FontSize', 12, 'FontWeight', 'bold');
