clear;
close all;
lambda = 3e-2; % Wavelength of the radar EM waves

%% Ground Truth Inputs

SNR_db = linspace(-30, 30, 1000); % SNR in db scale

SNR = 10.^(SNR_db./20); % SNR in linear scale




PRT = 1e-3;    % Pulse Repetition Time
% v_min = 2;     % Minimum Doppler velocity requierd in the spectrum
% v_max = 4;     % Maximum Doppler velocity required in the spectrum
mu = 3;        % Mean Doppler velocity for ground truth model
sigma = 0.2; % Standard deviation of velocity for ground truth model

beta_wind = eps;  % Azimuthal angle direction for the wind

BW_deg = 1;

BW = BW_deg*pi/180;
phi_0 = 0;
phi_end = 2*pi;

%% Use this when the radar has non-zero angular velocity
% Omega_rpm = linspace(1, 60, 120); % Angular velocity of radar in rpm
Omega_rpm = 1;
% Phi = phi_0:BW:phi_end;
% % 
% mean_Phi = mean([Phi(1:end-1); Phi(2:end)]);
mean_Phi = beta_wind;

%% Use this when the radar has 0 angular velocity
% Omega_rpm = 0;
% mean_Phi = beta_wind;
% N_burst = 256;
% X_ = rand(1, 1024);
% Theta_ = rand(1, 1024);

%% Signal generator and Doppler processing
for l = 1:length(SNR)
    for i = 1:length(Omega_rpm)
       Omega = Omega_rpm(i) .* 2.*pi./60; % Angular velocity of the radar beam in rad/sec
        for k = 1:length(mean_Phi)
                if Omega ~= 0
                    T = 2.*BW/Omega;
                else
                    T = N_burst .* PRT;
                end

                time_axis = eps:PRT:T; % Time axis in terms of multiples of PRT
                hits_scan(i) = length(time_axis); % length of time axis

                delta_v(i) = lambda/(2*hits_scan(i)*PRT);% velocity resolution
                v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
                vel_axis(i).axis = linspace(-v_amb,v_amb,hits_scan(i));% velocity axis

                Nifft = hits_scan(i); % Number of DTFT points for ifft


    %             beta = beta_wind - time_axis .* Omega; % Angle exerted by the radar beam with the wind direction at each time step
    %             noise(i).n = 1/SNR .* rand(1, hits_scan); % Noise based on the input SNR in linear scale
                beta = beta_wind - mean_Phi(k);



    %             X = X_(1:hits_scan(i));  % Random number generator for the simulator amplitude for all the velocities
    %             Theta = Theta_(1:hits_scan(i)) .* 2 .* pi; % Random phase generator for all the velocities (uncorelated with X in previous line)
    %             
                X(l).x = rand(1, hits_scan(i));  % Random number generator for the simulator amplitude for all the velocities
                Theta(l).theta = rand(1, hits_scan(i)) .* 2 .* pi; % Random phase generator for all the velocities (uncorelated with X in previous line)
            

                [Signal(i).sig(l, k, :)] = Doppler_spectrum_TD(vel_axis(i).axis, mu, sigma, Nifft, SNR(l), lambda, X(1).x, Theta(1).theta); 
                % Time domain signal generator without considering any rotation speed
                % of the radar beam


                Signal(i).sig_with_Omega(l, k, :) = abs(squeeze(Signal(i).sig(l, k, :))) .* exp(1j * angle(squeeze(Signal(i).sig(l, k, :))).* cos(beta));% + noise(i).n;
                % Use the rotation angle at each time step as an exponent (cos(beta))
                % to find the actual time domain signal for a rotating radar
    %             figure(101)
    %             hold on;
    % 
    %             txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
    %             plot(time_axis, abs(Signal(i).sig_with_Omega(k, :)), 'DisplayName',txt);
        end
    end
end

%% Processing of the simulated ground truth above 

for l = 1:length(SNR)
    for i = 1:length(Omega_rpm)
        for k = 1:length(mean_Phi)

                Nfft = hits_scan(i);

                Signal(i).doppler(l, k, :) = 1./sqrt(Nfft) .* (fftshift(fft(squeeze(Signal(i).sig_with_Omega(l, k, :))))); % Find the Doppler of the signal of the rotating radar



    %             [val, index] = max(squeeze(abs(Signal(i).doppler(k, :))));
    %             vel_max(i, k) = vel_axis(i).axis(index);

                PT_integrand = abs(squeeze(Signal(i).doppler(l, k, :))).^2 .* delta_v(i);
                PT = sum(PT_integrand); % Total power of the Doppler Spectrum

                v_mean_integrand = vel_axis(i).axis.' .* abs(squeeze(Signal(i).doppler(l, k, :))).^2 .* delta_v(i);
                v_mean(l, i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

                v_spread_integrand = (vel_axis(i).axis.' - squeeze(v_mean(l, i, k))).^2 .* abs(squeeze(Signal(i).doppler(l, k, :))).^2  .* delta_v(i);
                v_spread(l, i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width

                
                
                %% Mean Doppler Velocity




    %             figure(102)
    %             hold on;
    % 
    %             txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
    %             plot(vel_axis(i).axis, abs(Signal(i).doppler(k, :)), 'DisplayName',txt); % Plot Doppler spectrum with respect to the rotation speed of radar beam
    %         figure(101);
    % 
    %         txt_k = ['\Phi = ', num2str(mean_Phi(k).*180/pi), ' degree'];
    %         hold on; plot(vel_axis(i).axis, (abs(Signal(i).doppler(k, :))), 'DisplayName',txt_k);
        end
    %     xlabel('Velocity [m/s]', 'FontSize', 16);
    %     ylabel('Spectrum Power', 'FontSize', 16)
    %     title('Doppler Spectrum', 'FontSize', 16);
    %     legend show;
    %     grid on;
    %     figure(102);
    %     txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
    %    hold on; plot(mean_Phi * 180/pi, vel_max(i, :), 'DisplayName',txt);
    end
end


for l = 1:length(SNR)
    for i = 1:length(Omega_rpm)
        for k = 1:length(mean_Phi)
            v_meanE(l, i, k) = abs(squeeze(v_mean(l, i, k)) - mu)./squeeze(v_mean(end, i, k)) * 100;
            v_spreadE(l, i, k) = abs(squeeze(v_spread(l, i, k)) - sigma)./squeeze(v_spread(end, i, k)) * 100;
        end
    end
end

%% Plot erros wrt SNR

figure; plot(SNR_db, squeeze(v_meanE(:, 1, 1)), 'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2)
xlabel('SNR [dB]', 'FontSize', 16);
ylabel('Error in Doppler mean velocity [%]', 'FontSize', 16)
title('Error in Doppler mean velocity', 'FontSize', 16);
% legend show;
grid on;


figure; plot(SNR_db, squeeze(v_spreadE(:, 1, 1)), 'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2)
xlabel('SNR [dB]', 'FontSize', 16);
ylabel('Error in Doppler spectrum width [%]', 'FontSize', 16)
title('Error in Doppler spectrum width', 'FontSize', 16);
% legend show;
grid on;

%% Plot random Doppler spectrum for any angle of azimuth and any rpm listed on top
% 
% rpm_index = 1;
% phi_index = 1;
% 
% figure(2);hold on;
% plot(vel_axis(rpm_index).axis, sqrt(abs(Signal(rpm_index).doppler(phi_index, :)))); grid on;

%% Plot of mean velocity and velocity spectrum width with respect to rotation speed and azimuth angles

% figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_mean.'); shading flat; colormap('jet'); colorbar;
% 
% ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
% zlabel('Mean Doppler velocity [m/s]', 'FontSize', 16);
% 
% 
% 
% figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_spread.'); shading flat; colormap('jet'); colorbar;
% 
% ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
% zlabel('Doppler velocity width [m/s]', 'FontSize', 16);
% 
% %% Plot the Errors in mean and spectrum width of Doppler
% 
% figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_meanE.'); shading flat; colormap('jet'); colorbar;
% 
% ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
% zlabel('Mean Doppler velocity error [%]', 'FontSize', 16);
% % zlim([0 100]);
% 
% figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_spreadE.'); shading flat; colormap('jet'); colorbar;
% 
% ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
% zlabel('Doppler velocity width error [%]', 'FontSize', 16);
% % zlim([0 100]);