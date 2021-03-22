clear;
close all;
lambda = 3e-2; % Wavelength of the radar EM waves

%% Input

SNR_db = 40; % SNR in db scale

SNR = 10^(SNR_db./20); % SNR in linear scale



% Omega_rpm = linspace(1, 60, 10);
Omega_rpm = 1e-3;
% Omega_rpm = 0; % Angular velocity of radar in rpm
PRT = 1e-3;    % Pulse Repetition Time
v_min = 2;     % Minimum Doppler velocity requierd in the spectrum
v_max = 6;     % Maximum Doppler velocity required in the spectrum

beta_wind = eps;  % Azimuthal angle direction for the wind

BW_deg = 0.5;

BW = BW_deg*pi/180;
phi_0 = eps;
phi_end = 2*pi;

Phi = phi_0:BW:phi_end;

mean_Phi_axis = mean([Phi(1:end-1); Phi(2:end)]);
mean_Phi = mean_Phi_axis(1);

%% Signal generator and Doppler processing
for i = 1:1000
   Omega = Omega_rpm .* 2.*pi./60; % Angular velocity of the radar beam in rad/sec
    
            T = 2.*BW/Omega;

            time_axis = eps:PRT:T; % Time axis in terms of multiples of PRT
            hits_scan = length(time_axis); % length of time axis

            delta_v(i) = lambda/(2*hits_scan*PRT);% velocity resolution
            v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
            vel_axis = linspace(-v_amb,v_amb,hits_scan);% velocity axis

            Nifft = hits_scan; % Number of DTFT points for ifft


            beta = beta_wind - time_axis .* Omega; % Angle exerted by the radar beam with the wind direction at each time step
            noise(i).n = 1/SNR .* rand(1, hits_scan); % Noise based on the input SNR in linear scale

            X = rand(1, hits_scan);  % Random number generator for the simulator amplitude for all the velocities
            Theta = rand(1, hits_scan) .* 2 .* pi; % Random phase generator for all the velocities (uncorelated with X in previous line)

            [Signal(i, :), Signal_f(i, :)] = Doppler_spectrum_TD(vel_axis, v_min, v_max, Nifft, SNR, lambda, X, Theta); 
            % Time domain signal generator without considering any rotation speed
            % of the radar beam


            Signal_with_Omega(i, :) = abs((Signal(i, :))) .* (exp(1j * angle(Signal(i, :)))).^(cos(beta));% + noise(i).n;
            % Use the rotation angle at each time step as an exponent (cos(beta))
            % to find the actual time domain signal for a rotating radar

            Doppler(i, :) = fft((Signal_with_Omega(i, :))); % Find the Doppler of the signal of the rotating radar
            
   
end

Doppler_mean = mean(Signal_f, 1);

figure; plot(vel_axis, abs(Doppler_mean));

xlabel('velocity [m/s]', 'FontSize', 16);
ylabel('Doppler Spectrum', 'FontSize', 16)
title('Doppler Spectrum', 'FontSize', 16);
legend show;
grid on;

