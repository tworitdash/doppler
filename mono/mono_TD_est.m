clear;
close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1.8;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = 6;

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

beta_wind_deg = 0;

beta_wind = beta_wind_deg .* pi/180;

mu = 3;

[s] = TD_generator(mu, lambda, beta_wind, phi_0, Omega, t, 0.2);


figure; plot(t, unwrap(angle(s)));

%%

sec_in_one_rot = round((2*pi)/BW);
N_one_rot = sec_in_one_rot .* hs;

s_rot_reshape = reshape(s, N_one_rot, n_rot);
% s_rot_reshape_meas = reshape(s_meas, N_one_rot, n_rot);


s_sectors = zeros(sec_in_one_rot, n_rot*hs);
% s_sectors_meas = zeros(sec_in_one_rot, n_rot*hs);

for k = 1:sec_in_one_rot
   for i = 1:n_rot
      s_sectors(k, (i - 1)*hs+1:i*hs) = s_rot_reshape((k - 1)*hs+1:k*hs, i);
%       s_sectors_meas(k, (i - 1)*hs+1:i*hs) = s_rot_reshape_meas((k - 1)*hs+1:k*hs, i);
   end
end

%% Sector wise time domain estimaiton of Doppler velocity


for k = 1:sec_in_one_rot
    s_sectors_der(k, :) = diff(unwrap(angle(s_sectors(k, :))))./PRT .* lambda/(4 * pi);
end 

phi_axis_one_rot = linspace(0, 2*pi, sec_in_one_rot);


figure; surface(1:hs*n_rot-1, phi_axis_one_rot*180/pi, s_sectors_der); shading flat; colormap('jet'); colorbar;
