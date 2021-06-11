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
BW_deg = 1.8;
PRT = 1e-3;
lambda = 3e-2;
v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity


%% Radar rotation and Azimuth angle axes

% Omega_rpm = linspace(1, 60, 4); % RPM axis for the rotation of the radar

Omega_rpm = 1; % use this when only one rotation speed is needed 

Phi_0_deg = 0; % start of the azimuth angle of observation in degree 
Phi_end_deg = 360; % End of the azimuth angle of observation in degree

BW = BW_deg*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle in radians 
phi_end = Phi_end_deg * pi/180; % End of azimuth angle in radians

Phi = phi_0:BW:phi_end; % Azimuth angle axis separated with Beamwidths

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell


%% Data creation and processing


            
     Omega = Omega_rpm .* 2 * pi ./ 60; % Angular speed in [radians/sec]
     
     
     T = BW/Omega; % Dwell time (Time for the radar beam to stay at one azimuth angle cell (beamwidth))

     time_axis = eps:PRT:T; % Time axis for one beamwidth 

     hits_scan = length(time_axis); % length of time axis for one beam width (also known as hits per scan)
            
     N = hits_scan .* length(mean_Phi); % When creating HD spectrum and signal, we need number of points = hits_scan \times (Number of beam widths in an entire rotation)
            
     if mod(N, 2) == 0 % If it is even, make it odd
         N = N + 1;
     end
     
     phi_axis = linspace(phi_0, phi_end, N); % phi_axis = phi_axis(1:end-1);

     time_axis_full = linspace(eps, PRT.*N, N); % High resolution time axis
     axis_full = linspace(-(N - 1)/2, (N - 1)/2 - 1, N)/(N - 1); 
     vel_axis_full = 2.*v_amb.*axis_full; % High definition velocity axis 
            
     [data, data_f_orig, ~, ~] = ...
              DS_simulatorV3(SNR, m0, mu, sigma, N, v_amb); % Generation of HD signal in time domain
          
     data_f = 1./sqrt(N) .* fftshift(fft(data)); % Freq domain of HD
     
     
      data_man = abs(data) .* exp(1j .* (angle(data)) .* cos(beta_wind - phi_axis));
      
      data_man_f  = 1./sqrt(N) .* fftshift(fft(data_man(1:N))); % Freq domain of manipulated
      
      
      
      
      
      data_man_re = abs(data_man) .* exp(1j .* (angle(data_man)) ./ (cos(phi_axis)));
      
      
      
      data_man_re_f = 1./sqrt(N) .* fftshift(fft(data_man_re));
      
      
      
%       figure; plot(vel_axis_full, db(abs(data_f))); hold on;  
      
      plot(vel_axis_full, db(abs(data_f_orig))); 
      
      hold on; plot(vel_axis_full, db(abs(data_man_f))); 
      
      hold on; plot(vel_axis_full, db(abs(data_man_re_f))); 
       
      legend({'HD spectrum', 'Manipulated with cosine Spectrum', 'Compensated from Manipulated spectrum'});  grid on;
        
        
%             
%       figure; plot(angle(data)); hold on; plot(angle(data_man_re))



%% Make sectors and process




data_man_sec = abs(data) .* exp(1j .* unwrap(angle(data)) .* cos(beta_wind - phi_axis));


for k = 1:length(mean_Phi)
    
    sig_doppler(k, :) = 1./sqrt(hits_scan) .* fftshift(fft(data_man_sec((k - 1) * hits_scan+1:k*hits_scan)));
    
    
end



vel_axis = linspace(-v_amb, v_amb, hits_scan);

figure; imagesc(vel_axis, mean_Phi.*180/pi, db(abs(sig_doppler)));


