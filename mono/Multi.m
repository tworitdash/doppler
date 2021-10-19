clear;
close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1.8;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
% Omega_rpm = linspace(1, 60, 4);
Omega_rpm = 1;
% Omega_rpm = 1;

lambda = 0.03;
%% 

for i = 1:length(Omega_rpm)
   
    
BW = BW_deg * pi/180;
Omega(i) = Omega_rpm(i) * 2*pi/60; % rotation speed in rad/s
Td = BW/Omega(i);
hs = round(Td/PRT);

hs_(i) = hs;

phi_0 = phi_0_deg * pi/180;

sec = round((n_rot*2*pi)/BW);

N = sec * hs;

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


mu = normrnd(mu_mean, 1, [1 hs]);
% mu = normrnd(4, 1, [1 N]);
n_sig = sqrt(0.01);

% mu = 3;


if mod(hs, 2) == 0
    vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
%      vel_axis_hs = linspace(-v_amb, v_amb, hs+1);
else
    vel_axis_hs = linspace(-v_amb, v_amb, hs);
end

for k = 1:length(phi_axis-1)
    t = 0:PRT:(hs - 1)*PRT;
    phi = linspace(phi_axis_1(k), phi_axis_1(k+1), hs);
    [s, SNR] = TD_generator(mu, lambda, beta_wind, phi_axis_1(k), Omega(i), t, n_sig);
    
    s_fft(k, :) = 1./sqrt(hs) .* fftshift(fft(s)); 
    
end


end

figure; surface(vel_axis_hs,phi_axis*180/pi, db(abs(s_fft).^2)/2); shading flat; 