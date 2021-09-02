clear;
close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1.8;
n_rot = 3;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = 1;

lambda = 0.03;
%% 

BW = BW_deg * pi/180;
Omega = Omega_rpm * 2*pi/60; % rotation speed in rad/s
Td = BW/Omega;
hs = round(Td/PRT);
phi_0 = phi_0_deg * pi/180;

sec = round((n_rot*2*pi)/BW);
N = sec * hs;
t = 0:PRT:(N - 1)*PRT;
phi_axis = phi_0:BW:phi_0+Omega*t(end);
% phi_axis = zeros(size(phi_axis));


%% variables for signal model

beta_wind_deg = 90;

beta_wind = beta_wind_deg .* pi/180;


mu = normrnd(3, 1, [1 N]);
n_sig = 0.1;

% mu = 3;
[s] = TD_generator(mu, lambda, beta_wind, phi_0, Omega, t, n_sig);

%% variables for real life scenario

% beta_wind_meas_deg = 0;
% 
% beta_wind_meas = beta_wind_meas_deg .* pi/180;
% 
% mu_meas = 3;
% 
% [s_meas] = TD_meas(mu_meas, lambda, beta_wind_meas, phi_0, Omega, t);


%% Time domain correction

% figure; plot(t(1:hs), real(s(1:hs)), 'LineWidth', 2); hold on; plot(t(1:hs), real(s_meas(1:hs)), '*');
% 
% figure; plot(t(1:hs), imag(s(1:hs)), 'LineWidth', 2); hold on; plot(t(1:hs), imag(s_meas(1:hs)), '*');
% 

%% Store the time domain data sector-wise

sec_in_one_rot = round((2*pi)/BW);
N_one_rot = sec_in_one_rot .* hs;

s_rot_reshape = reshape(s, N_one_rot, n_rot);
% s_rot_reshape_meas = reshape(s_meas, N_one_rot, n_rot);


s_sectors = zeros(sec_in_one_rot, n_rot*hs);
s_sectors_meas = zeros(sec_in_one_rot, n_rot*hs);

for k = 1:sec_in_one_rot
   for i = 1:n_rot
      s_sectors(k, (i - 1)*hs+1:i*hs) = s_rot_reshape((k - 1)*hs+1:k*hs, i);
%       s_sectors_meas(k, (i - 1)*hs+1:i*hs) = s_rot_reshape_meas((k - 1)*hs+1:k*hs, i);
   end
end

phi_axis_one_rot = linspace(0, 2*pi, sec_in_one_rot);

% figure; surface(1:hs*n_rot, phi_axis_one_rot*180/pi, real(s_sectors)); shading flat; colormap('jet');
% 
% figure; surface(1:hs*n_rot, phi_axis_one_rot*180/pi, real(s_sectors_meas)); shading flat; colormap('jet');

%%


N = length(t);

v_amb = lambda/(4 * PRT);

if mod(hs, 2) == 0
    vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
else
    vel_axis_hs = linspace(-v_amb, v_amb, hs);
end

[S1, F1, Ti1, P1] = spectrogram(s, hs, 0, hs , N*1/60);


S1 = S1./max(S1(:));
S1_norm = fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));



figure;
surface(vel_axis_hs, phi_axis*180/pi, S1_norm_db); shading flat;
colormap('jet');
colorbar;

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title('Doppler Spectrum');


%% Calculation of mean
diff_v = diff(vel_axis_hs); 
dv = diff_v(1); % velocity resolution

for k = 1:size(S1_norm, 1)
   
    PT = sum(abs(S1_norm(k, :)).^2 .* dv);
    v_mean_i(k, :) = vel_axis_hs .* abs(S1_norm(k, :)).^2 .* dv;
    v_mean(k) = 1./PT .* sum(v_mean_i(k, :));
    
    v_spread_i(k, :) = (vel_axis_hs - v_mean(k)).^2 .* abs(S1_norm(k, :)).^2 .* dv;
    v_spread(k) = sqrt(1./PT .* sum(v_spread_i(k, :)));
    
end
% mu_der = (mu(2) - mu(1))./(t(2) - t(1));
% beta_wind_der = beta_wind(2) - beta_wind(1);
% v_max = (-1/Omega) .* (mu_der .* sin(- phi_0) + mu .* (beta_wind_der - Omega).*cos( - phi_0));
% txt = ['Angle integration = ', num2str(N_BW*BW_deg), ' [deg]'];
% figure(101); hold on; plot(phi * 180/pi, v_mean, 'LineWidth', 2, 'DisplayName', txt); grid on;
figure(101); hold on; plot(phi_axis * 180/pi, v_mean, 'LineWidth', 2); grid on; % hold on; plot(phi_0 + Omega .* t, v_max); 
legend;
title('Mean Doppler velocity')
xlabel('Angle [deg]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('v_{mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');

figure(102); hold on; plot(phi_axis * 180/pi, sqrt(v_spread), 'LineWidth', 2); grid on; % hold on; plot(phi_0 + Omega .* t, v_max); 
legend;
title('Doppler spectrum width')
xlabel('Angle [deg]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('v_{mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');

%% LS method to find all three component of velocity

theta = 0; % elevation

T = 2 * pi * n_rot/Omega;



u = 1./cos(theta) .* 2/(Omega .* T) .* sum(v_mean .* cos(phi_axis) .* BW);

v = 1./cos(theta) .* 2/(Omega .* T) .* sum(v_mean .* sin(phi_axis) .* BW);

U = sqrt(abs(u).^2 + abs(v).^2);
beta_wind_re_deg = angle(u + 1j .* v) .* 180/pi;

