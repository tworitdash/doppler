%% Input to the simulator 
clear;
close all;

Omega_rpm = linspace(1, 60, 5);
% Omega_rpm = 0;
SNR_db = inf;
BW_deg = 1; 

Phi_0_deg = 0;
Phi_end_deg = 360;

beta_wind = pi/2;
PRT = 1e-3;
lambda = 3e-2;

mu = 5;
sigma = 0.2;


BW = BW_deg*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle
phi_end = Phi_end_deg * pi/180; % End of azimuth angle

Phi = phi_0:BW:phi_end;

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell


N_iteration = 1;

for i = 1:length(Omega_rpm)
%     if Omega_rpm == 0
%         mean_Phi = beta_wind;
%     end
    
    for k = 1:length(mean_Phi)
    
        for l = 1:N_iteration % Monte Carlo Iterations

            [Signal(i).sig(l, k, :), Signal(i).sig_doppler(l, k, :), Signal(i).sig_with_Omega(l, k, :), hits_scan(i), delta_v(i), vel_axis(i).axis, time_axis(i).axis] ...
                = Simulator_with_rot(Omega_rpm(i), BW, SNR_db, mean_Phi(k), beta_wind, PRT, lambda, mu, sigma);

            Signal(i).doppler(l, k, :) = 1./sqrt(hits_scan(i)) .* (fftshift(fft(squeeze(Signal(i).sig_with_Omega(l, k, :))))); % ./cos(beta_wind - mean_Phi(k)); % Find the Doppler of the signal of the rotating radar
            Signal(i).doppler_norm(l, k, :) =  abs(squeeze(Signal(i).doppler(l, k, :)))./ max(abs(squeeze(Signal(i).doppler(l, k, :))));
            PT_integrand = abs(squeeze(Signal(i).doppler(l, k, :))).^2 .* delta_v(i);
            PT = sum(PT_integrand); % Total power of the Doppler Spectrum

            v_mean_integrand = vel_axis(i).axis.' .* abs(squeeze(Signal(i).doppler(l, k, :))).^2 .* delta_v(i);
            v_mean(i, l, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

            v_spread_integrand = (vel_axis(i).axis.' - v_mean(i, l, k)).^2 .* abs(squeeze(Signal(i).doppler(l, k, :))).^2 .* delta_v(i);
            v_spread(i, l, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
            
        end
        v_mean_error(i, k) = sqrt(sum((squeeze(v_mean(i, :, k)) - mu).^2)/N_iteration);
        v_mean_ik(i, k) = mean(v_mean(i, :, k));
        Signal(i).doppler_norm_avg(k, :) = mean(squeeze(Signal(i).doppler_norm(:, k, :)), 1);
    end

    figure; surface(vel_axis(i).axis, mean_Phi * 180/pi, (Signal(i).doppler_norm_avg)); shading flat; colorbar; %colormap('jet'); 
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    title(['Normalized Doppler spectrum at ', num2str(Omega_rpm(i)), 'RPM'], 'FontSize', 16);
    
    
    figure; surface(vel_axis(i).axis, mean_Phi * 180/pi, (abs(squeeze(Signal(i).sig(1, :, :))))); shading flat; colorbar; %colormap('jet'); 
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    title(['Normalized Doppler spectrum at ', num2str(Omega_rpm(i)), 'RPM'], 'FontSize', 16);


end


% figure; surface(Omega_rpm, mean_Phi, v_mean_error.'); shading flat; colorbar;colormap('jet') % Plot the error in mean Doppler
% 
% ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
% zlabel('Error in Mean Doppler velocity [m/s]', 'FontSize', 16);
% 
% figure; surface(Omega_rpm, mean_Phi, v_mean_ik.'); shading flat; colorbar; colormap('jet'); % Plot the mean Doppler
% 
% ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
% zlabel('Mean Doppler velocity [m/s]', 'FontSize', 16);

