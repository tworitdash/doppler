
clear;
close all;
lambda = 3e-2; % Wavelength of the radar EM waves

%% Input

SNR_db = 30; % SNR in db scale

SNR = 10^(SNR_db./20); % SNR in linear scale



% Omega_rpm = linspace(0, 60, 2);
Omega_rpm = 0; % Angular velocity of radar in rpm
PRT = 1e-3;    % Pulse Repetition Time
v_min = 2;     % Minimum Doppler velocity requierd in the spectrum
v_max = 6;     % Maximum Doppler velocity required in the spectrum

beta_wind = eps;  % Azimuthal angle direction for the wind

%% Signal generator and Doppler processing
for k = 1:1000
for i = 1:length(Omega_rpm)
   

    Omega = Omega_rpm(i) .* 2.*pi./60; % Angular velocity of the radar beam in rad/sec
    
    if Omega ~= 0
        T = 2.*pi/Omega/8;  % Time period for one full rotation of the radar beam
    else
        T = 2^12.*PRT;    % Total time of observation when the radar is static. It has to be given manually as it can't be determined from the rotation speed
    end
    
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
    
    [Signal(i).sig, Vol(k, :)] = Doppler_spectrum_TD(vel_axis, v_min, v_max, Nifft, SNR, lambda, X, Theta); 
    % Time domain signal generator without considering any rotation speed
    % of the radar beam
    
    
    Signal(i).sig_with_Omega = abs((Signal(i).sig)) .* (exp(1j * angle(Signal(i).sig))).^(cos(beta));% + noise(i).n;
    % Use the rotation angle at each time step as an exponent (cos(beta))
    % to find the actual time domain signal for a rotating radar
    
    Signal(i).doppler = fft((Signal(i).sig_with_Omega)); % Find the Doppler of the signal of the rotating radar
    
%     figure(101)
%     hold on;
    
%     txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
%     plot(vel_axis, abs(Signal(i).doppler), 'DisplayName',txt); % Plot Doppler spectrum with respect to the rotation speed of radar beam
    
end
end

figure; plot(vel_axis, mean(Vol, 1));
% xlabel('Velocity [m/s]');
% title('Doppler Vel');
% legend show;
% grid on;
% 
% 
% figure;
% plot(Omega_rpm, delta_v, 'LineWidth', 2); % Plot of Doppler velocity resolution with respect to the rotation speed of the radar beam.
% xlabel('Rotation speed [rpm]');
% ylabel('Doppler velocity resolution [m/s]');
% title('Doppler velocity resolution vs rotation speed of radar in rpm')
% % legend show;
% grid on;


%% Pinsky Algorithm 


data = Signal(1).sig; % Time domain radar complex signal


a = zeros(hits_scan, 1); % Initialization of a (The parameter to be estimated)
Theta_model = ones(hits_scan,1); % Initialization  of the variance estimation of a 
data_model = zeros(hits_scan, 1); % Initialization of the estimate of the radar data retrieval

gamma = zeros(hits_scan, 1); % Initialization of Gamma which helps in reduction in the analytical model of the estimator 

gamma(1) = abs(data(1)).^(-2); % Filling the first entry of gamma as 1/square_root(first entry in measurement data)
P = ones(hits_scan, 1); % Covariance matrix initialization for the estimated data (data_model)
Q = zeros(hits_scan, 1); % Initialization of the dynamic SNR 

% 
% data_model(1, :) = data(1, :);

tic;

for  k = 2:hits_scan
    %disp(k)
    
    gamma(k) = gamma(k - 1) ./ (1 + gamma(k - 1) .* abs(data_model(k - 1)).^2); % update for Gamma
    a(k) = a(k - 1) + gamma(k) .* conj(data_model(k -1)) .* (data(k) - a(k - 1) .* data_model(k - 1)); % update for a
    
    
    
%     Q(k - 1) = sum(abs(data(1:k-1))).^2 ./ (mn.^2 .* sum(1:k - 1)) - 1;
    Q(k - 1) = sum(abs(data(1:k-1))).^2 ./ (sum(abs(noise(1).n(1:k - 1)))).^2 - 1; % SNR till all measurements before time 'k'
    P(k) = (Q(k - 1) .* (1 - abs(a(k)).^2) + P(k - 1) .* abs(a(k)).^2) ./ ...
        (1 + Q(k - 1) .* (1 - abs(a(k)).^2) + P(k - 1) .* abs(a(k)).^2); % Update for P
    data_model(k) = a(k) .* data_model(k - 1) + P(k) .* (data(k) - a(k) .* data_model(k - 1)); % Update for estimated data 
    
end
db(Q)
omega = unwrap(angle(a)); % Mean Doppler velocity


var_omega_i = zeros(10000, hits_scan); % Initialization of the Doppler spectrum width
for l = 1:10000
    var_omega_i(l, :) = 4 .* (-1).^(l + 1) .* (abs(a)).^l;  
end

var_omega = pi^2/3 - sum(var_omega_i, 1); % Doppler spectrum width

v = omega .* lambda ./ 2;
time_adap = toc;


figure; 
hold on;
plot(1:hits_scan, v, '*'); hold on; plot(1:hits_scan, var_omega);
ylim([0 5]);

xlabel('Hits scan')
title('velocity of the targets [m/s]')
