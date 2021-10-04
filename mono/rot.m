clear;
close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1.8;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = linspace(1, 60, 4);
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

% sec_in_one_rot = round((2*pi)/BW);
% N_one_rot = sec_in_one_rot .* hs;
% 
% s_rot_reshape = reshape(s, N_one_rot, n_rot);
% % s_rot_reshape_meas = reshape(s_meas, N_one_rot, n_rot);
% 
% 
% s_sectors = zeros(sec_in_one_rot, n_rot*hs);
% s_sectors_meas = zeros(sec_in_one_rot, n_rot*hs);
% 
% for k = 1:sec_in_one_rot
%    for i = 1:n_rot
%       s_sectors(k, (i - 1)*hs+1:i*hs) = s_rot_reshape((k - 1)*hs+1:k*hs, i);
% %       s_sectors_meas(k, (i - 1)*hs+1:i*hs) = s_rot_reshape_meas((k - 1)*hs+1:k*hs, i);
%    end
% end
% 
% phi_axis_one_rot = linspace(0, 2*pi, sec_in_one_rot);

% figure; surface(1:hs*n_rot, phi_axis_one_rot*180/pi, real(s_sectors)); shading flat; colormap('jet');
% 
% figure; surface(1:hs*n_rot, phi_axis_one_rot*180/pi, real(s_sectors_meas)); shading flat; colormap('jet');

%%


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
% 

%% Calculation of mean
diff_v = diff(vel_axis_hs); 
dv = diff_v(1); % velocity resolution

if i > 1
    clear v_mean_i;
    clear v_spread_i;
end


for k = 1:size(S1_norm, 1)
   
    PT = sum(abs(S1_norm(k, :)).^2 .* dv);
    v_mean_i(k, :) = vel_axis_hs .* abs(S1_norm(k, :)).^2 .* dv;
    v_mean(i, l, k) = 1./PT .* sum(v_mean_i(k, :));
    
    v_spread_i(k, :) = (vel_axis_hs - v_mean(i, l, k)).^2 .* abs(S1_norm(k, :)).^2 .* dv;
    v_spread(i, l, k) = sqrt(1./PT .* sum(v_spread_i(k, :)));
    
end
end
end
% mu_der = (mu(2) - mu(1))./(t(2) - t(1));
% beta_wind_der = beta_wind(2) - beta_wind(1);
% v_max = (-1/Omega) .* (mu_der .* sin(- phi_0) + mu .* (beta_wind_der - Omega).*cos( - phi_0));
% txt = ['Angle integration = ', num2str(N_BW*BW_deg), ' [deg]'];
% figure(101); hold on; plot(phi * 180/pi, v_mean, 'LineWidth', 2, 'DisplayName', txt); grid on;
%% 
idx_phi = 1;

txt = ['SNR = ', num2str(db(SNR)/2), ' dB', ', BeamWidth = ', num2str(BW_deg), '^{\circ}', ' \phi = ', num2str(phi_axis(idx_phi) .* 180/pi), '^{\circ}'];
% 
% figure;  imagesc(Omega_rpm, mu_mean, squeeze(v_mean(:, :, idx_phi)).'); colormap('jet'); colorbar; shading flat;% , '-', 'LineWidth', 2, 'DisplayName', txt); grid on; % hold on; plot(phi_0 + Omega .* t, v_max); 
% % legend;
% title(['Mean Doppler velocity', txt])
% xlabel('\Omega [rpm]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('\mu_{True mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
% 
% xlim([Omega_rpm(1) Omega_rpm(end)]);
% ylim([mu_mean(1) mu_mean(end)]);
% 
% figure;  imagesc(Omega_rpm, mu_mean, squeeze(v_spread(:, :, idx_phi)).'); colormap('jet'); colorbar; shading flat; %, '-', 'LineWidth', 2, 'DisplayName', txt); % hold on; plot(phi_0 + Omega .* t, v_max); 
% % legend;
% title(['Doppler spectrum width', txt])
% xlabel('\Omega [rpm]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('\mu_{True mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
% 


for i = 1:length(Omega_rpm)
    name = ['\Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];
    figure(101);
    hold on; plot(mu_mean, squeeze(v_mean(i, :, 1)), 'DisplayName', name);
    figure(102);
    hold on; plot(mu_mean, squeeze(v_spread(i, :, 1)), 'DisplayName', name)
end
figure(101)
grid on;
lgd = legend;
lgd.FontSize = 14;
lgd.FontWeight = 'bold';
title(['Mean Doppler velocity', txt], 'FontSize', 12, 'FontWeight', 'bold')
ylabel('\mu_{re} [m.s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('\mu_{True mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');

xlim([mu_mean(1) mu_mean(end)]);
ylim([mu_mean(1) mu_mean(end)]);


figure(102)
grid on;
lgd = legend;
lgd.FontSize = 14;
lgd.FontWeight = 'bold';
title(['Doppler spectrum width', txt], 'FontSize', 12, 'FontWeight', 'bold')
ylabel('\sigma_{re} [m.s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('\mu_{True mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');

%% 
idx_mean = 50;

txt = ['SNR = ', num2str(db(SNR)/2), ' dB', ', BeamWidth = ', num2str(BW_deg), '^{\circ}', ' \mu_{mean} = ', num2str((mu_mean(idx_mean))), ' [m.s^{-1}]'];


figure; imagesc(Omega_rpm, phi_axis .* 180/pi, squeeze(v_mean(:, idx_mean, :)).'); colormap('jet'); colorbar; shading flat;% , '-', 'LineWidth', 2, 'DisplayName', txt); grid on; % hold on; plot(phi_0 + Omega .* t, v_max); 
% legend;
title(['Mean Doppler velocity', txt])
xlabel('\Omega [rpm]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\phi [\circ]', 'FontSize', 12, 'FontWeight', 'bold');

xlim([Omega_rpm(1) Omega_rpm(end)]);
ylim([phi_axis(1)*180/pi phi_axis(end)*180/pi]);

figure; imagesc(Omega_rpm, phi_axis .* 180/pi, squeeze(v_spread(:, idx_mean, :)).'); colormap('jet'); colorbar; shading flat; %, '-', 'LineWidth', 2, 'DisplayName', txt); % hold on; plot(phi_0 + Omega .* t, v_max); 
% legend;
title(['Doppler spectrum width', txt])
xlabel('\Omega [rpm]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\phi [\circ]', 'FontSize', 12, 'FontWeight', 'bold');

xlim([Omega_rpm(1) Omega_rpm(end)]);
ylim([phi_axis(1)*180/pi phi_axis(end)*180/pi]);

%%

%% 
idx_omega = 1;

txt = ['SNR = ', num2str(db(SNR)/2), ' dB', ', BeamWidth = ', num2str(BW_deg), '^{\circ}', ', \Omega = ', num2str((Omega_rpm(idx_omega))), ' [rpm]'];


figure; surface(mu_mean, phi_axis .* 180/pi, squeeze(v_mean(idx_omega, :, :)).'); colormap('jet'); colorbar; shading flat;% , '-', 'LineWidth', 2, 'DisplayName', txt); grid on; % hold on; plot(phi_0 + Omega .* t, v_max); 
% legend;
title(['Mean Doppler velocity', txt])
xlabel('\mu_{True mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\phi [\circ]', 'FontSize', 12, 'FontWeight', 'bold');

xlim([mu_mean(1) mu_mean(end)]);
ylim([phi_axis(1)*180/pi phi_axis(end)*180/pi]);

figure; surface(mu_mean, phi_axis .* 180/pi, squeeze(v_spread(idx_omega, :, :)).'); colormap('jet'); colorbar; shading flat; %, '-', 'LineWidth', 2, 'DisplayName', txt); % hold on; plot(phi_0 + Omega .* t, v_max); 
% legend;
title(['Doppler spectrum width', txt])
xlabel('\mu_{True mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\phi [\circ]', 'FontSize', 12, 'FontWeight', 'bold');

xlim([mu_mean(1) mu_mean(end)]);
ylim([phi_axis(1)*180/pi phi_axis(end)*180/pi]);


%% VAD technique 


