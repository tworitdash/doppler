clc; 
clear; 
% close all;
c0 = 3e8; % Speed of EM waves
% Radar specification 

lambda = 3e-2; % Wavelength 

vel = 5; % velocity of wind 

beta_wind = eps; % direction of wind in terms of azimuth

%% SNR



SNR_db = Inf; % SNR In db scale
% SNR_db = linspace(-30, 30, 10);

SNR = 10.^(SNR_db/20); % SNR in linear scale



%% Use when radar rotation speed is non zero

% Omega_rpm = linspace(1, 60, 3);
Omega_rpm = 1;
BW_deg = 1;

BW = BW_deg*pi/180;
phi_0 = eps;
phi_end = 2*pi;
Phi = phi_0:BW:phi_end;

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell
%% Use when radar rotation speed is 0
% 
% % Omega_rpm = 0;
% mean_Phi = beta_wind;
% N_burst = 256;

%% radar PRT
PRT = 1e-3; % Pulse repetition time

figure;
% for l = 1:length(SNR)
    for i = 1:length(Omega_rpm)

       Omega = Omega_rpm(i) .* 2*pi/60;
        for k = 1:length(mean_Phi)

                if Omega ~= 0
                    T = BW/Omega;
                else
                    T = N_burst .* PRT;
                end

                time_axis = eps:PRT:T; % Time axis in terms of multiples of PRT
                hits_scan(i) = length(time_axis); % length of time axis

                delta_v(i) = lambda/(2*hits_scan(i)*PRT);% velocity resolution
                v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
                vel_axis(i).axis = linspace(-v_amb,v_amb,hits_scan(i));% velocity axis

                Nifft = hits_scan(i); % Number of DTFT points for ifft


                Vel_with_Omega_0(k) = vel .* cos(beta_wind - mean_Phi(k)); % True value of radial velocity looked from an angle mean_Phi(k)


                Signal(i).sig(k, :) = data_simulator_BW(time_axis, mean_Phi(k), beta_wind, SNR, lambda, vel); % Signal generator output

                Signal(i).doppler(k, :) = fftshift(fft(squeeze(Signal(i).sig(k, :)), hits_scan(i))); % FFT is done with same number of points as the length of the signal
                
                PT_integrand = abs(Signal(i).doppler(k, :)).^2 .* delta_v(i);
                PT = sum(PT_integrand); % Total power of the Doppler Spectrum
            
                v_mean_integrand = vel_axis(i).axis .* abs(Signal(i).doppler(k, :)).^2 .* delta_v(i);
                v_mean(i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 
            
                v_spread_integrand = (vel_axis(i).axis - v_mean(i, k)).^2 .* abs(Signal(i).doppler(k, :)).^2 .* delta_v(i);
                v_spread(i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
            
           
            if k == 52
                figure(16); hold on; plot(time_axis, angle(Signal(i).sig(k, :)))
            end
                
%                 [val, index] = max(squeeze(Signal(i).doppler(l, k, :)));
%                 vel_max(l, i, k) = vel_axis(index);
%                 
%                 e_v(l, i, k) = abs(squeeze(vel_max(l, i, k) - Vel_with_Omega_0(k))./Vel_with_Omega_0(k)) .* 100;
%                 
%                 figure(101);
%                 if k == 1
%                     hold on; plot(vel_axis, abs(squeeze(Signal.doppler(l, k, :))));
%                 end
    %             txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
    %             plot(vel_axis, abs(Signal(i).doppler), 'DisplayName',txt); % Plot of the Doppler Power spectrum w.r.t the velocity axis

        end
        
    figure; surface(vel_axis(i).axis, mean_Phi * 180/pi, abs(Signal(i).doppler)); shading flat; colorbar; %colormap('jet');
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    title(['Normalized Doppler spectrum at ', num2str(Omega_rpm(i)), 'RPM'], 'FontSize', 16);
    figure(101)
    
    txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
    hold on; plot(mean_Phi * 180/pi, v_mean(i, :), '*', 'DisplayName',txt);
    end
%     figure(102)
%     txt = ['SNR = ', num2str(SNR_db(l)), ' dB' ];
%     hold on; plot(mean_Phi * 180/pi, squeeze(v_mean(l, i, :)), 'DisplayName',txt);
% end

xlabel('Azimuth Angle \Phi [deg]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Mean Doppler Velocity observed [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
title('Maximum Doppler velocity observed', 'FontSize', 16, 'FontWeight', 'bold');
h = legend;
set(h,'FontSize',16, 'FontWeight', 'bold', 'Location','north');
grid on;

% figure(103);
% surface(SNR_db, mean_Phi .* 180/pi, squeeze(e_v(:, 1, :)).'); shading flat; colormap('jet'); colorbar;
% 
% figure(104);
% plot(SNR_db, squeeze(e_v(:, 1, 1)).');

% figure;
% % plot(Omega_rpm, delta_v, 'LineWidth', 2); 
% xlabel('Rotation speed [rpm]');
% ylabel('Doppler velocity resolution [m/s]');
% title('Doppler velocity resolution vs rotation speed of radar in rpm')
% % legend show;
% grid on;


%% Plot spectrum at any angle and any rpm
% figure;
% Omega_rpm_i = [1 2 3 4 5 6 7];
% mean_Phi_i = [180];
% 
% for i = 1:length(Omega_rpm_i)
%     for k = 1:length(mean_Phi_i)
%         txt = ['\Omega = ', num2str(Omega_rpm(Omega_rpm_i(i))), ' rpm'];
%         hold on; 
%         plot(vel_axis(i).axis, abs(Signal(i).doppler(mean_Phi_i(k), :)));
%     end
% end
% 
% xlabel('velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Spectrum', 'FontSize', 16, 'FontWeight', 'bold')
% title('Velocity spectrum', 'FontSize', 16, 'FontWeight', 'bold');
% h = legend;
% set(h,'FontSize',16, 'FontWeight', 'bold', 'Location','north');
% grid on;

%% Plot at 51.5 degree

% figure(301); hold on; plot(vel_axis(1).axis, abs(Signal(1).doppler(52, :))/max(abs(Signal(1).doppler(52, :))))

figure(301); hold on; plot(vel_axis(1).axis, abs(Signal(1).doppler(52, :)))