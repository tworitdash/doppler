clear;
close all;
lambda = 3e-2; % Wavelength of the radar EM waves

%% Ground Truth Inputs

SNR_db = -8; % SNR in db scale

SNR = 10^(SNR_db./20); % SNR in linear scale




PRT = 1e-3;    % Pulse Repetition Time
% v_min = 2;     % Minimum Doppler velocity requierd in the spectrum
% v_max = 4;     % Maximum Doppler velocity required in the spectrum
mu = 3;        % Mean Doppler velocity for ground truth model
sigma = 0.1; % Standard deviation of velocity for ground truth model

beta_wind = eps;  % Azimuthal angle direction for the wind

BW_deg = 1;

BW = BW_deg*pi/180;
phi_0 = 0;
phi_end = 2*pi;

%% Use this when the radar has non-zero angular velocity
% Omega_rpm = linspace(0.2, 60, 120); % Angular velocity of radar in rpm
Omega_rpm = 0.2;
Phi = phi_0:BW:phi_end;
% 
mean_Phi = mean([Phi(1:end-1); Phi(2:end)]);
% mean_Phi = pi/2;

%% Use this when the radar has 0 angular velocity
% Omega_rpm = 0;
% mean_Phi = beta_wind;
% N_burst = 256;
X_ = rand(1, 1024);
Theta_ = rand(1, 1024);

%% Signal generator and Doppler processing
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
            noise(i).n = 1/SNR .* rand(1, hits_scan(1)); % Noise based on the input SNR in linear scale
            beta = beta_wind - mean_Phi(k);
            
          
            
%             X = X_(1:hits_scan(i));  % Random number generator for the simulator amplitude for all the velocities
%             Theta = Theta_(1:hits_scan(i)) .* 2 .* pi; % Random phase generator for all the velocities (uncorelated with X in previous line)
%             
            X1(i).x = rand(1, hits_scan(i));  % Random number generator for the simulator amplitude for all the velocities
            Theta1(i).theta = rand(1, hits_scan(i)) .* 2 .* pi; % Random phase generator for all the velocities (uncorelated with X in previous line)
            
%             idx = 1:hits_scan(1);
%             idxq = linspace(min(idx), max(idx), hits_scan(i));
%             X = interp1(idx, X1(1).x, idxq, 'linear');
%             Theta = interp1(idx, Theta1(1).theta, idxq, 'linear');
%          
            
            [Signal(i).sig(k, :), Signal(i).sig_doppler(k, :)] = Doppler_spectrum_TD(vel_axis(i).axis, mu, sigma, Nifft, SNR, lambda, X1(i).x, Theta1(i).theta); 
            % Time domain signal generator without considering any rotation speed
            % of the radar beam

            PT_integrand_simulator = abs(Signal(i).sig_doppler(k, :)).^2 .* delta_v(i);
            PT_simulator = sum(PT_integrand_simulator); % Total power of the Doppler Spectrum
            
            v_mean_integrand_simulator = vel_axis(i).axis .* abs(Signal(i).sig_doppler(k, :)).^2 .* delta_v(i);
            v_mean_simulator(i, k) = sum(v_mean_integrand_simulator) ./ PT_simulator; % Mean Doppler velocity 
            
            v_spread_integrand_simulator = (vel_axis(i).axis - v_mean_simulator(i, k)).^2 .* abs(Signal(i).sig_doppler(k, :)).^2 .* delta_v(i);
            v_spread_simulator(i, k) = sqrt(sum(v_spread_integrand_simulator)./ PT_simulator); % Doppler spectrum width
            
            Signal(i).sig_with_Omega(k, :) = abs((Signal(i).sig(k, :))) .* exp(1j * angle(Signal(i).sig(k, :)).* cos(beta));% + noise(i).n;
            % Use the rotation angle at each time step as an exponent (cos(beta))
            % to find the actual time domain signal for a rotating radar
%             figure(101)
%             hold on;
% 
%             txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
%             plot(time_axis, abs(Signal(i).sig_with_Omega(k, :)), 'DisplayName',txt);
    end
end

%% Processing of the simulated ground truth above 


for i = 1:length(Omega_rpm)
    for k = 1:length(mean_Phi)
            
            Nfft = hits_scan(i);
            
            Signal(i).doppler(k, :) = 1./sqrt(Nfft) .* (fftshift(fft((Signal(i).sig_with_Omega(k, :))))); % Find the Doppler of the signal of the rotating radar
            
%             Signal(i).doppler(k, :) = (fftshift(fft((Signal(i).sig_with_Omega(k, :))))); % Find the Doppler of the signal of the rotating radar
            
%             [val, index] = max(squeeze(abs(Signal(i).doppler(k, :))));
%             vel_max(i, k) = vel_axis(i).axis(index);
            
            PT_integrand = abs(Signal(i).doppler(k, :)).^2 .* delta_v(i);
            PT = sum(PT_integrand); % Total power of the Doppler Spectrum
            
            v_mean_integrand = vel_axis(i).axis .* abs(Signal(i).doppler(k, :)).^2 .* delta_v(i);
            v_mean(i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 
            
            v_spread_integrand = (vel_axis(i).axis - v_mean(i, k)).^2 .* abs(Signal(i).doppler(k, :)).^2 .* delta_v(i);
            v_spread(i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
            
            v_meanE(i, k) = abs(v_mean(i, k) - v_mean_simulator(i, k))./v_mean_simulator(i, k) * 100;
            v_spreadE(i, k) = abs(v_spread(i, k) - v_spread_simulator(i, k))./v_spread_simulator(i, k) * 100;
            
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


% xlabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
% ylabel('Maximum Velocity observed [m/s]', 'FontSize', 16)
% title('Maximum Doppler velocity observed', 'FontSize', 16);
% legend show;
% grid on;


%% Plot random Doppler spectrum for any angle of azimuth and any rpm listed on top
% 
% rpm_index = 1;
% phi_index = 1;
% 
% figure(2);hold on;
% plot(vel_axis(rpm_index).axis, sqrt(abs(Signal(rpm_index).doppler(phi_index, :)))); grid on;

%% Plot of mean velocity and velocity spectrum width with respect to rotation speed and azimuth angles

figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_mean.'); shading flat; colormap('jet'); colorbar;

ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
zlabel('Mean Doppler velocity [m/s]', 'FontSize', 16);



figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_spread.'); shading flat; colormap('jet'); colorbar;

ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
zlabel('Doppler velocity width [m/s]', 'FontSize', 16);

%% Plot the Errors in mean and spectrum width of Doppler

figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_meanE.'); shading flat; colormap('jet'); colorbar;

ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
zlabel('Mean Doppler velocity error [%]', 'FontSize', 16);
% zlim([0 100]);

figure; surface(Omega_rpm, mean_Phi .* 180/pi, v_spreadE.'); shading flat; colormap('jet'); colorbar;

ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
zlabel('Doppler velocity width error [%]', 'FontSize', 16);
% zlim([0 100]);


%% Jordi algorithm 

data = Signal(1).sig(1, :);

a = zeros(1, hits_scan(1));
Theta_model = ones(1, hits_scan(1));
data_model = zeros(1, hits_scan(1));
% data_model(1) = data(1);
data_model_var = ones(1, hits_scan(1));
gamma = zeros(1, hits_scan(1));

gamma(1) = abs(data(1)).^(-2);
P = zeros(1, hits_scan(1));
P(1)= 1;
Q = zeros(1, hits_scan(1));

% 
% data_model(1, :) = data(1, :);

tic;

for  k = 1:hits_scan(1)-1
    %disp(k)
    Q(k) = sum(abs(data(1:k))).^2 ./ (sum(noise(1).n(1:k)).^2 - 1);
    
    P_ = Q(k) .* (1 - abs(a(k)).^2) + P(k) .* (abs(a(k)).^2);
    P(k + 1) = P_./(1 + P_);
 
    data_model(k + 1) = a(k) .* data_model(k) + P(k + 1) .* (data(k + 1) - a(k) .* data_model(k));
    
    gamma(k + 1) = gamma(k) ./ (1 + gamma(k) .* abs(data_model(k)).^2);
    a(k + 1) = a(k) + gamma(k + 1) .* conj(data_model(k)) .* (data(k + 1) - a(k) .* data_model(k));
end



omega = (angle(a));

% F = pn.^2./(1 - a .* exp(-1j .* omega));
v = omega .* lambda / 2;
time_adap = toc;


figure; plot(1:hits_scan(1), v, 'LineWidth', 2); figure; plot(1:hits_scan(1), abs(a), 'LineWidth', 2);
figure; plot(1:hits_scan(1), abs(data), 'LineWidth', 2); hold on; plot(1:hits_scan(1), abs(data_model), 'LineWidth', 2);
figure; plot(1:hits_scan, db(Q));
% v_max = 10;
% v_abs = abs(v);
% v_abs(v_abs>v_max) = v_max;
% 
% data_model_abs = abs(data_model);
% data_model_max = 70;
% data_model_abs(data_model_abs>data_model_max) = data_model_max;
% 
% figure; surface(range*1e-3, 1:hits_scan, abs(data)); shading flat; colorbar, colormap('jet'); 
% 
% xlabel('Range [km]')
% ylabel('Hits scan')
% title('Raw data in time domain')
% 
% figure; surface(range*1e-3, 1:hits_scan, data_model_abs); shading flat; colorbar, colormap('jet');caxis([0 70]);
% 
% xlabel('Range [km]')
% ylabel('Hits scan')
% title('Data after filtering, Estimate of the compelx data')
% 
% figure; surface(range*1e-3, 1:hits_scan, v_abs); shading flat; colorbar, colormap('jet'); caxis([1 10]);
% 
% xlabel('Range [km]')
% ylabel('Hits scan')
% title('velocity of the targets [m/s]')
% 
% 
% figure;
% imagesc(range*1e-3,[],abs(data)); colorbar;
% xlabel('Range [km]')
% ylabel('Hits per scan')
% title('Raw data 2D plot')
% 
% %6. Make Doppler
% tic;
% data_Doppler = fftshift(fft(data,[],1), 1);
% time_fft = toc;
% figure;
% imagesc(vel,range*1e-3,abs(data_Doppler).'); shading flat; colorbar; colormap('jet');
% % surface(vel,range*1e-3,abs(data_Doppler).'); shading flat; colorbar; colormap('jet');
% xlabel('Velocity [m/s]')
% ylabel('Range [km]')
% title('Range Doppler Map')