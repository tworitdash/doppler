
%% RADAR constants and wind direction

beta_wind = eps;

SNR_db = Inf;
SNR = 10^(SNR_db/20);

PRT = 1e-3;
lambda = 3e-2;
n = 10;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity

[data, data_f] = DS_simulator(SNR, mu, sigma, n, v_amb); % One time simulation of the spectrum and the signal with fine sampling

%% 

BW_deg = 1;

Phi_0_deg = 0;
Phi_end_deg = 360;

BW = BW_deg*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle
phi_end = Phi_end_deg * pi/180; % End of azimuth angle

Phi = phi_0:BW:phi_end;

mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell

Omega_rpm = linspace(1, 60, 5); % RPM axis for the rotation of the radar

%% resampling based on the rotation of the radar

for i = 1:length(Omega_rpm)
    
    Omega = Omega_rpm(i) .* 2 * pi ./ 60;
    
    T = BW/Omega;
    
    time_axis(i).axis = eps:PRT:T;
    
    hits_scan_(i) = length(time_axis(i).axis); % length of time axis
    hits_scan(i) = 2^(nextpow2(hits_scan_(i)) - 1); % hits scan for Doppler processing
    vel_axis(i).axis = linspace(-v_amb, v_amb, hits_scan(i));
    
    for k = 1:length(mean_Phi)
        beta(k) = beta_wind - mean_Phi(k);
        Signal(i).sig(k, :) = abs(data) .* exp(1j .* unwrap(angle(data)) .* cos(beta(k)));
        Signal(i).doppler(k, :) = 1./sqrt(hits_scan(i)) .* fftshift(fft(Signal(i).sig(k, :)));
    end
    
    figure; surface(vel_axis(i).axis, mean_Phi * 180/pi, abs(Signal(i).doppler)); shading flat; colorbar; %colormap('jet'); 
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    title(['Normalized Doppler spectrum at ', num2str(Omega_rpm(i)), 'RPM'], 'FontSize', 16);
    
 end