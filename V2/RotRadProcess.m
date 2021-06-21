clear; 

%% Ground Truth
beta_wind = eps; % wind direction
m0 = 1;
mu = 4;
sigma = 1;

%% Radar parameters 

PRT = 1e-3; 
BW_deg = 1.8;
SNR_db = 30;
lambda = 0.03;
Omega_rpm = 1;
Phi_0_deg = 0;
Phi_end_deg = 360;

%% Signal simulator output
[sig, data, data_f, N, hs, phi_axis, mean_Phi, vel_axis_full, time_axis_full, v_amb] = RotSimulator(PRT, BW_deg, SNR_db, lambda, m0, mu, sigma, Omega_rpm, Phi_0_deg, Phi_end_deg, beta_wind);

figure; plot(time_axis_full, abs(ifftshift(data)));
figure; plot(time_axis_full, real(ifftshift(data))); hold on; plot(time_axis_full, imag(ifftshift(data)));


figure; plot(vel_axis_full, db(abs(fftshift(fft(sig)))));
figure; plot(vel_axis_full, db(abs(fftshift(fft(data)))));
figure; plot(vel_axis_full, db(abs(data_f)));
%% Processing

%% Spectogram

N_BW = 1; % Number of beamwidths to consider for Doppler processing
M = round(length(mean_Phi)/N_BW);
phi = linspace(phi_axis(1), phi_axis(end), M); % Angle of the sectors

dv = lambda/(2*hs*PRT);

% if mod(hs, 2) == 0
%     hs = hs + 1;
% end

vel_axis_hs = linspace(-v_amb, v_amb, hs);

[S1, F1, Ti1, P1] = spectrogram(sig, hs*N_BW, 0, hs, N*Omega_rpm/(60));


S1 = S1./max(S1(:));
S1_norm = fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));

F1_doppler = F1 - N/2;

V_doppler = lambda .* F1_doppler ./ 2;

figure;
imagesc(vel_axis_hs, phi * 180/pi, S1_norm_db); shading flat;
colormap('jet');
colorbar;

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title('Doppler Spectrum');

%% 1D plots of velocity with azimuth at different rotation speeds but a given Beamwidth
 

for k = 1:length(phi) % For each direction in azimuth 
          
      PT_integrand = abs(squeeze(abs(S1_norm(k, :)))).^2 .* dv;
      PT = sum(PT_integrand); % Total power of the Doppler Spectrum

      v_mean_integrand = vel_axis_hs .* abs(squeeze(S1_norm(k, :))).^2 .* dv;
      v_mean(k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

      v_spread_integrand = (vel_axis_hs - v_mean(k)).^2 .* abs(squeeze(S1_norm(k, :))).^2 .* dv;
      v_spread(k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
                

end
            
figure;

txt = ['\Omega = ', num2str(Omega_rpm), ' [rpm]'];
hold on; plot(phi * 180/pi, squeeze(v_mean), 'DisplayName', txt);

h = legend;
grid on;
ylabel('Mean Doppler velocity [m/s]', 'FontSize', 12);
xlabel('Azimuth \Phi [deg]', 'FontSize', 12)
set(h,'FontSize',12, 'FontWeight', 'bold', 'Location','north');
title(['Mean Doppler velocity when ', ' SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 12);

figure;

txt = ['\Omega = ', num2str(Omega_rpm), ' [rpm]'];
hold on; plot(phi * 180/pi, squeeze(v_spread), 'DisplayName', txt);

   
h = legend;
grid on;
ylabel('Doppler spectrum width [m/s]', 'FontSize', 12);
xlabel('Azimuth \Phi [deg]', 'FontSize', 12)
set(h,'FontSize',12, 'FontWeight', 'bold', 'Location','north');
title(['Doppler spectrum width when ', ' SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 12);

