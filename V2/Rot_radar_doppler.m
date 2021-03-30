clear;
close all;

%% Configuration 0 is off, 1 in on

% Want plot for Error vs SNR? 
SNR_enable = 1;

% Want plot for Erorr vs Omega and Phi?
OP_enable_error = 0;

% Want Doppler Azimuth plot?

DS_Azimuth_plots = 0; % Doppler Azimuth Plot

% Want to run Monte Carlo simulation?
MC_enable = 1;

% Want to plot mean Doppler and Doppler spectrum width?
OP_enable = 0;

%%

if MC_enable == 1
    n_MC = 32; 
else
    n_MC = 1;
end


if SNR_enable == 1

    SNR_i_db = -30; 
    SNR_f_db = 30; 
    n_SNR = 10; % Number of points wanted in the SNR axis
    SNR_db = linspace(SNR_i_db, SNR_f_db, n_SNR);
    SNR = 10.^(SNR_db/20);

else
    n_SNR = 1;
    SNR_db = Inf;
    SNR = 10^(SNR_db/20);
end


%% RADAR constants and wind direction

beta_wind = pi/2; % wind direction
mu = 5;
sigma = 0.2;



PRT = 1e-3;
lambda = 3e-2;
n = 2^10;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity

% for l = 1:n_MC
%     for s = 1:n_SNR
%         [data_(l, s, :), data_f_(l, s, :)] = DS_simulator(SNR(s), mu, sigma, n, v_amb); % One time simulation of the spectrum and the signal with fine sampling
%     end
% end

% data = squeeze(mean(data_, 1));
% data_f = squeeze(mean(data_f_, 1));

% n_ = 4;
% vel_ = linspace(-v_amb, v_amb, n);
% vel__ = linspace(-v_amb, v_amb, n_);
% 
% figure; plot(vel_, abs(data_f)); hold on; plot(vel__, 1/sqrt(n) .* abs(fftshift(fft(data, n_))));

%% 

BW_deg = 1;

Phi_0_deg = 0;
Phi_end_deg = 360;

BW = BW_deg*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle
phi_end = Phi_end_deg * pi/180; % End of azimuth angle

Phi = phi_0:BW:phi_end;

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell

Omega_rpm = linspace(1, 12, 4); % RPM axis for the rotation of the radar
% Omega_rpm = 160;
% Omega_rpm = 60;

%% resampling based on the rotation of the radar
for l = 1:n_MC % Monte Carlo iterations
    
    for s = 1:n_SNR
        
        [data(l, s, :), data_f(l, s, :)] = DS_simulator(SNR(s), mu, sigma, n, v_amb);
        
        for i = 1:length(Omega_rpm)

            Omega = Omega_rpm(i) .* 2 * pi ./ 60;

            T = BW/Omega;

            time_axis(i).axis = eps:PRT:T;

            hits_scan_(i) = length(time_axis(i).axis); % length of time axis
            hits_scan(i) = 2^(nextpow2(hits_scan_(i)) - 1); % hits scan for Doppler processing
            vel_axis(i).axis = linspace(-v_amb, v_amb, hits_scan(i));
            delta_v(i) = lambda/(2*hits_scan(i)*PRT);

            for k = 1:length(mean_Phi)
               beta(k) = beta_wind - mean_Phi(k);
    %             beta(k) = eps;
                Signal(i).sig(l, s, k, :) = abs(squeeze(data(l, s, :))) .* exp(1j .* unwrap(angle(squeeze(data(l, s, :)))) .* cos(beta(k)));
                Signal(i).doppler(l, s, k, :) = 1./sqrt(hits_scan(i)) .* fftshift(fft(Signal(i).sig(l, s, k, :), hits_scan(i)));

                PT_integrand = abs(squeeze(Signal(i).doppler(l, s, k, :))).^2 .* delta_v(i);
                PT = sum(PT_integrand); % Total power of the Doppler Spectrum

                v_mean_integrand = vel_axis(i).axis.' .* abs(squeeze(Signal(i).doppler(l, s, k, :))).^2 .* delta_v(i);
                v_mean_l(l, s, i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

                v_spread_integrand = (vel_axis(i).axis.' - v_mean_l(l, s, i, k)).^2 .* abs(squeeze(Signal(i).doppler(l, s, k, :))).^2 .* delta_v(i);
                v_spread_l(l, s, i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
            end
            
    %         figure; surface(vel_axis(i).axis, mean_Phi * 180/pi, abs(squeeze(Signal(i).doppler(s, :, :)))); shading flat; colorbar; %colormap('jet'); 
    %         ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    %         xlabel('Velocity axis [m/s]', 'FontSize', 16)
    %         zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    %         title(['Normalized Doppler spectrum at ', num2str(Omega_rpm(i)), 'RPM', ', SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 16);
        end

    end
end

% v_mean = squeeze(mean(v_mean_l, 1));
% v_spread = squeeze(mean(v_spread_l, 1));

%% Error calculations


for s = 1:n_SNR
   for i = 1:length(Omega_rpm)
      for k = 1:length(mean_Phi)
           v_mean(s, i, k) = sum(v_mean_l(:, s, i, k))/n_MC;
           v_spread(s, i, k) = sum(v_spread_l(:, s, i, k))/n_MC;
           v_mean_e(s, i, k) = sqrt(sum((v_mean_l(:, s, i, k) - mu).^2)/n_MC);
           v_spread_e(s, i, k) = sqrt(sum((v_spread_l(:, s, i, k) - sigma).^2)/n_MC);
      end
   end
end


%% Plot Doppler spectrum (Azimuth vs Doppler Velocity)


if DS_Azimuth_plots == 1
    SI = 1; % Index for SNR
    OI = 1;   % Index for Omega
    Plot2DDoppler(vel_axis(OI).axis, mean_Phi, Signal, SI, OI, SNR_db, Omega_rpm);
end
%% Plot Erros

if SNR_enable == 1
    OI = 1; % Index for Omega
    PI = 1; % Index for Phi
    PlotDopplerSNR(SNR_db, v_mean_e, v_spread_e, OI, PI, Omega_rpm, mean_Phi);
end

%% 2D plots of mean Doppler and Doppler spread and Erros

if OP_enable == 1 || OP_enable_error == 1
    SI = 1; % Index of the SNR axis
    PlotDopplerOP(OP_enable, OP_enable_error, SNR_db, SI, mean_Phi, Omega_rpm, v_mean, v_spread, v_mean_e, v_spread_e);
end

%% 1D plots of velocity with azimuth at different rotation speeds but a given Beamwidth


figure;
SI = 1; % Index of the SNR axis
for i = 1:length(Omega_rpm)
   txt = ['\Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];
   hold on; plot(mean_Phi * 180/pi, squeeze(v_mean(SI, i, :)), 'DisplayName', txt);
end
h = legend;
grid on;
ylabel('Mean Doppler velocity [m/s]', 'FontSize', 16);
xlabel('Azimuth \Phi [deg]', 'FontSize', 16)
set(h,'FontSize',16, 'FontWeight', 'bold', 'Location','north');
title(['Mean Doppler velocity when ', ' BW = ', num2str(BW_deg), ' [deg] and SNR = ', num2str(SNR_db(SI)), ' dB'], 'FontSize', 16);


