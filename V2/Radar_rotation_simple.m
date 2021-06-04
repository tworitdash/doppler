clear;
close all;

%% Ground Truth 


beta_wind = eps; % wind direction
m0 = 1; % Amplitiude
mu = 5; % Mean velocity
sigma = 0.2; % Standard deviation of wind velocity

%% Radar variables 


SNR_db = 30;
SNR = 10^(SNR_db/10);
BW_deg = 1;
PRT = 1e-3;
lambda = 3e-2;
v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity


%% Radar rotation and Azimuth angle axes

Omega_rpm = linspace(1, 60, 4); % RPM axis for the rotation of the radar

% Omega_rpm = 40.33; % use this when only one rotation speed is needed 

Phi_0_deg = 0; % start of the azimuth angle of observation in degree 
Phi_end_deg = 360; % End of the azimuth angle of observation in degree

BW = BW_deg*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle in radians 
phi_end = Phi_end_deg * pi/180; % End of azimuth angle in radians

Phi = phi_0:BW:phi_end; % Azimuth angle axis separated with Beamwidths

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell


%% Data creation and processing


for i = 1:length(Omega_rpm) % Loop over each angular velocity of the radar
            
     Omega = Omega_rpm(i) .* 2 * pi ./ 60; % Angular speed in [radians/sec]

     T = BW/Omega; % Dwell time (Time for the radar beam to stay at one azimuth angle cell (beamwidth))

     time_axis = eps:PRT:T; % Time axis for one beamwidth 
            

     hits_scan = length(time_axis); % length of time axis for one beam width (also known as hits per scan)
            
     N(i) = hits_scan .* length(Phi); % When creating HD spectrum and signal, we need number of points = hits_scan \times (Number of beam widths in an entire rotation)
            
     if mod(N(i), 2) ~= 0 % If it is even, make it odd
         N(i) = N(i) + 1;
     end
            
            
%% This if clause is only to decide the number of points required for Doppler processing for every azimuthal angle
     if hits_scan < 32 % (Make it 0 if no zero padding is required)
         hits_scan_1(i) = 64 + 1; % For Doppler processing when hits per scan is too small [zero padding purposes]
     else
         hits_scan_1(i) = hits_scan; 
     end
%% 
     time_axis_full = linspace(eps, PRT.*N(i), N(i)); % High resolution time axis
     axis_full = linspace(-(N(i) - 1)/2, (N(i) - 1)/2 - 1, N(i))/(N(i) - 1); 
     vel_axis_full = 2.*v_amb.*axis_full; % High definition velocity axis 
            
     [data(i).data, data(i).data_f, ~, ~] = ...
              DS_simulatorV3(SNR, m0, mu, sigma, N(i), v_amb); % Generation of HD signal in time domain
          
     % (i -> for each rotation speed, m -> For different BW, l ->
     % Monte Carlo axis, s -> SNR axis, k -> For each azimuth angle) 
         
     vel_axis = linspace(-v_amb, v_amb, hits_scan_1(i)); % velocity axis to show the Doppler spectrum and calculate the moments
            
     delta_v = lambda/(2*hits_scan_1(i)*PRT); % velocity resolution 

     for k = 1:length(mean_Phi) % For each direction in azimuth 
                
         beta_scan = beta_wind - (Phi(k) + Omega .* time_axis);
            
                 
         Signal(i).sig(k, :) = hamming(hits_scan).' .* (abs(squeeze(data(i).data((k - 1)*hits_scan + 1: k*hits_scan)))...
                .* exp(1j .* unwrap(angle(squeeze(data(i).data((k - 1)*hits_scan + 1: k*hits_scan)))) .* cos(beta_scan)));
             
         Signal(i).doppler(k, :) = 1./sqrt(hits_scan_1(i)) .* fftshift(fft(Signal(i).sig(k, :), hits_scan_1(i))); % To be used for moments calculation
                
                                              
%        figure; plot(vel_axis, db(abs(squeeze(Signal(i).doppler(k, :)).^2))/2, '-o');               
%        hold on;  plot(vel_axis_full, db(abs(squeeze(data(i).data_f)).^2)/2, '-.'); grid on;
 
                
              
           %% Moments calculation
           
          PT_integrand = abs(squeeze(Signal(i).doppler(k, :))).^2 .* delta_v;
          PT = sum(PT_integrand); % Total power of the Doppler Spectrum

          v_mean_integrand = vel_axis .* abs(squeeze(Signal(i).doppler(k, :))).^2 .* delta_v;
          v_mean(i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

          v_spread_integrand = (vel_axis - v_mean(i, k)).^2 .* abs(squeeze(Signal(i).doppler(k, :))).^2 .* delta_v;
          v_spread(i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
                

     end
            
end

%% Plot Doppler spectrum (Azimuth vs Doppler Velocity)
OI = 1; % Index for Omega

Plot2DDopplerV2(hits_scan_1(OI), mean_Phi, v_amb, Signal, OI, SNR_db, Omega_rpm);


%% 1D plots of velocity with azimuth at different rotation speeds but a given Beamwidth


figure;


for i = 1:length(Omega_rpm)
   txt = ['\Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];
   hold on; plot(mean_Phi * 180/pi, squeeze(v_mean(i, :)), 'DisplayName', txt);
end
h = legend;
grid on;
ylabel('Mean Doppler velocity [m/s]', 'FontSize', 12);
xlabel('Azimuth \Phi [deg]', 'FontSize', 12)
set(h,'FontSize',12, 'FontWeight', 'bold', 'Location','north');
title(['Mean Doppler velocity when ', ' SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 12);
figure;
for i = 1:length(Omega_rpm)
   txt = ['\Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];
   hold on; plot(mean_Phi * 180/pi, squeeze(v_spread(i, :)), 'DisplayName', txt);
end
h = legend;
grid on;
ylabel('Doppler spectrum width [m/s]', 'FontSize', 12);
xlabel('Azimuth \Phi [deg]', 'FontSize', 12)
set(h,'FontSize',12, 'FontWeight', 'bold', 'Location','north');
title(['Doppler spectrum width when ', ' SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 12);

