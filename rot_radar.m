clc; 
clear; 
% close all;
c0 = 3e8; % Speed of EM waves
% Radar specification 

lambda = 3e-2; % Wavelength 

vel = 2; % velocity of wind 

beta_wind = eps; % direction of wind in terms of azimuth

SNR_db = Inf; % SNR In db scale
SNR = 10^(SNR_db/20); % SNR in linear scale

Omega_rpm = linspace(0, 60, 3);
% Omega_rpm = 0;
PRT = 1e-3; % Pulse repetition time


figure;
for i = 1:length(Omega_rpm)
   

    Omega = Omega_rpm(i) .* 2.*pi./60; % Rotation speed in radian per second
    
    if Omega ~= 0
        T = 2*pi/Omega/10; % Time period for one rotation 
    else
        T = 256.*PRT; % 256 bursts are used in the case of a stationary radar
    end
    
    time_axis = eps:PRT:T; % Time axis for one rotation
    
    hits_scan = length(time_axis); % Number of bursts 
    
    delta_v(i) = lambda/(2*hits_scan*PRT);% velocity resolution
    v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
    vel_axis = linspace(-v_amb,v_amb,hits_scan);% velocity axis


    [Signal(i).sig] = data_simulator(time_axis, Omega, beta_wind, SNR, lambda, vel); % Signal generator output
    
    Signal(i).doppler = fftshift(fft((Signal(i).sig), hits_scan)); % FFT is done with same number of points as the length of the signal
    
    hold on;
    
    txt = ['\Omega = ', num2str(Omega_rpm(i)), ' rpm' ];
    plot(vel_axis, abs(Signal(i).doppler), 'DisplayName',txt); % Plot of the Doppler Power spectrum w.r.t the velocity axis
    
end

xlabel('Velocity [m/s]');
title('Doppler retrieval');
legend show;
grid on;

figure;
% plot(Omega_rpm, delta_v, 'LineWidth', 2); 
xlabel('Rotation speed [rpm]');
ylabel('Doppler velocity resolution [m/s]');
title('Doppler velocity resolution vs rotation speed of radar in rpm')
% legend show;
grid on;

%% 
% figure; 
% surf(Omega_rpm, vel_axis, abs(Signal.doppler)); shading flat; colormap('jet'); colorbar;
% xlabel('Velocity [m/s]');
% xlabel('\Omega [rpm]');
% title('Doppler Vel');

% %% Pinsky Algorithm 
% 
% 
% data = Signal;
% 
% 
% a = zeros(hits_scan, 1);
% Theta_model = ones(hits_scan,1);
% data_model = zeros(hits_scan, 1);
% 
% gamma = zeros(hits_scan, 1);
% 
% gamma(1) = abs(data(1)).^(-2);
% P = ones(hits_scan, 1);
% Q = zeros(hits_scan, 1);
% 
% % 
% % data_model(1, :) = data(1, :);
% 
% tic;
% 
% for  k = 2:hits_scan
%     %disp(k)
%     
%     gamma(k) = gamma(k - 1) ./ (1 + gamma(k - 1) .* abs(data_model(k - 1)).^2);
%     a(k) = a(k - 1) + gamma(k) .* conj(data_model(k -1)) .* (data(k) - a(k - 1) .* data_model(k - 1));
%     
%     
%     
% %     Q(k - 1) = sum(abs(data(1:k-1))).^2 ./ (mn.^2 .* sum(1:k - 1)) - 1;
%     Q(k - 1) = sum(abs(data(1:k-1))).^2 ./ (sum(abs(noise(1:k - 1)))).^2 - 1;
%     P(k) = (Q(k - 1) .* (1 - abs(a(k)).^2) + P(k - 1) .* abs(a(k)).^2) ./ ...
%         (1 + Q(k - 1) .* (1 - abs(a(k)).^2) + P(k - 1) .* abs(a(k)).^2);
%     data_model(k) = a(k) .* data_model(k - 1) + P(k) .* (data(k) - a(k) .* data_model(k - 1));
%     
% end
% db(Q)
% omega = angle(a);
% 
% 
% var_omega_i = zeros(10000, hits_scan);
% for l = 1:10000
%     var_omega_i(l, :) = 4 .* (-1).^(l + 1) .* (abs(a)).^l; 
% end
% 
% var_omega = pi^2/3 - sum(var_omega_i, 1);
% 
% v = omega .* lambda ./ 2;
% time_adap = toc;
% 
% 
% 
% 
% figure; 
% hold on;
% plot(1:hits_scan, omega); hold on; plot(1:hits_scan, var_omega);
% ylim([0 5]);
% 
% xlabel('Hits scan')
% title('velocity of the targets [m/s]')
% 
% 
% %% Dynamic model


