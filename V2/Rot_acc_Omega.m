clear; 
close all;

%% Wind ground truth

beta_wind = 0;
theta = 0;
mu = 2;
sigma = 0.1;

%% Radar

PRT = 1e-3;

% Omega_rpm = linspace(20, 23, 4);
% Omega_rpm = 1;
Omega_rpm = ones(1, 10) .* 40;

BW_deg = 1;
BW = BW_deg * pi / 180;

Phi = eps:BW:2*pi;

lambda = 3e-2;
n = 2^10;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity

%% SNR
% SNR_db = linspace(-70, 70, 100);
SNR_db = 30;
SNR = 10.^(SNR_db/10);
n_MC = 16;

%% Loop for data accumulation 

for m = 1:n_MC
   for s = 1:length(SNR_db)
        for i = 1:length(Omega_rpm)
            Omega = Omega_rpm(i) .* (2 .* pi) ./ 60;
            Td = BW/Omega;
            time_axis_scan = eps:PRT:Td;
            
            hits_scan_time = length(time_axis_scan);
            
            hits_scan(i) = 2^(nextpow2(hits_scan_time) - 1);
            
            T_residual = (2.*pi - BW)/Omega;
            time_axis_residual = eps:PRT:T_residual;
            hits_scan_residual(i) = 2^(nextpow2(length(time_axis_residual)) - 1);
            sweep_total(i) = hits_scan(i) + hits_scan_residual(i);
           
            
            for k = 1:length(Phi)-1
                beta_scan_hd = beta_wind - linspace(Phi(k), Phi(k + 1), n);
                Signal(i).sig_m(m, s, k, :) = [DS_simulatorV3_with_az(SNR(s), 1, mu, sigma, n, v_amb, hits_scan(i), beta_scan_hd, theta) zeros(1, hits_scan_residual(i))];
            end
        end
   end
end

for i = 1:length(Omega_rpm)
    Signal(i).sig(:, :, :) = mean(squeeze(Signal(i).sig_m), 1);
end

All_scan = [0 cumsum(2.^(nextpow2(sweep_total) - 1))];

for l = 1:length(All_scan)
    for s = 1:length(SNR_db)
        for i = 1:length(Omega_rpm)
            for k = 1:length(Phi)-1
                sig(s, k, All_scan(i)+1:All_scan(i +1)) = squeeze(Signal(i).sig(s, k, 1:All_scan(2)));
            end
        end
    end
end

%% Plot time domain
SI = 1;
PI = 1;
time_axis = eps:PRT:(All_scan(end) - 1)* PRT;
figure; plot(time_axis.', real(squeeze(sig(SI, PI, :))), 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]); 
hold on; plot(time_axis.', imag(squeeze(sig(SI, PI, :))), 'LineWidth', 2, 'color', [0.0780, 0.6350, 0.1840]); grid on;

title(['Time domain signal at SNR = ', num2str(SNR_db), ' dB and \Phi = ' num2str((Phi(PI) + Phi(PI + 1))/2 .* 180/pi), ' [deg]']);
legend({'Real', 'Imaginary'})



%% Frequency domain analysis

N = 2^(nextpow2(All_scan(end)));
Axis = linspace(-N/2, N/2-1, N);
vel_axis = 2 .* v_amb .* Axis ./ N;


sig_doppler = 1./sqrt(N) .* fftshift(fft(sig, N, 3), 3);
%% Plot the Doppler spectrum
SI = 1;
figure; imagesc(vel_axis, Phi(1:length(Phi)-1).*180/pi, abs(squeeze(sig_doppler(SI, :, :))).^2); shading flat; colormap('jet'); colorbar; title(['Doppler Spectrum at SNR = ', num2str(SNR_db), ' dB'])
%% Plot for a specific angle 

PI = 1;
SI = 1;
figure; plot(vel_axis, abs(squeeze(sig_doppler(SI, PI, :))).^2, 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]); grid on;

%% Mean Doppler velocity 
 
delta_v = lambda/(2*N*PRT);

for s = 1:length(SNR_db)
    for k = 1:length(Phi)-1
        PT_integrand = abs(squeeze(sig_doppler(s, k, :))).^2 .* delta_v;
        PT = sum(PT_integrand); % Total power of the Doppler Spectrum

        v_mean_integrand = vel_axis.' .* abs(squeeze(sig_doppler(s, k, :))).^2 .* delta_v;
        v_mean(s, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

        v_spread_integrand = (vel_axis.' - v_mean(s, k)).^2 .* abs(squeeze(sig_doppler(s, k, :))).^2 .* delta_v;
        v_spread(s, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
    end
end

%% Plot mean 
SI = 1;
figure; plot(Phi(1:length(Phi) - 1) .* 180/pi, v_mean(SI, :), 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840] ); title('Mean Doppler Velocity', 'FontSize', 16); grid on;

% figure; plot(Phi(1:length(Phi) - 1) .* 180/pi, v_spread, 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]); title('Doppler Spectrum Width', 'FontSize', 16); grid on;
%% Finding the direction of wind flow

[val, idx] = max(v_mean(:, 1:181), [], 2);
phi_wind = Phi(idx) .* 180/pi;

figure; plot(SNR_db, phi_wind, 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]); grid on;