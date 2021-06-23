function [sig, data, data_f, N, hits_scan, phi_axis, mean_Phi, vel_axis_full, time_axis_full, v_amb] = RotSimulator(PRT, BW_deg, SNR_db, lambda, m0, mu, sigma, Omega_rpm, Phi_0_deg, Phi_end_deg, beta_wind)


%Rotating radar time domain data simulator. 
%  


%% Radar variables 

SNR = 10^(SNR_db/10);

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity


BW = BW_deg*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle in radians 
phi_end = Phi_end_deg * pi/180; % End of azimuth angle in radians

Phi = phi_0:BW:phi_end; % Azimuth angle axis separated with Beamwidths

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell


%% Data creation and processing


 % Loop over each angular velocity of the radar
            
 Omega = Omega_rpm .* 2 * pi ./ 60; % Angular speed in [radians/sec]

 T = BW/Omega; % Dwell time (Time for the radar beam to stay at one azimuth angle cell (beamwidth))

 time_axis = 0:PRT:T; % Time axis for one beamwidth 
            

 hits_scan = length(time_axis); % length of time axis for one beam width (also known as hits per scan)
 if mod(hits_scan, 2) == 0 % If it is even, make it odd
     hits_scan = hits_scan + 1;
 end
                
 N = hits_scan .* length(mean_Phi); % When creating HD spectrum and signal, we need number of points = hits_scan \times (Number of beam widths in an entire rotation)
            
 if mod(N, 2) == 0 % If it is even, make it odd
     N = N + 1;
 end
     
 phi_axis = linspace(phi_0, phi_end, N); % phi_axis = phi_axis(1:end-1);
%  mu = mu .* cos(phi_axis);
%  sigma = sigma .* cos(phi_axis);
            
 time_axis_full = linspace(0, PRT.*N, N); % High resolution time axis
vel_axis_full = linspace(-v_amb, v_amb, N);
            
 [data, data_f, ~, ~] = ...
              DS_simulatorV3(SNR, m0, mu, sigma, N, v_amb); % Generation of HD signal in time domain
 sig = abs(data) .* exp(1j .* unwrap(angle(data)) .* sin(eps + beta_wind - phi_axis) ./ (eps - Omega .* time_axis_full));
          
         

            