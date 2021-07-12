clear; 
close all;


lambda = 0.03;



Omega_rpm = 1; % in RPM
Omega = 2*pi/60 * Omega_rpm; % In rad/s
BW_deg = 1.8; % beam width in degree
BW = BW_deg * pi/180; % beam width in radian

v_amb = 7.5;


PRT = 1e-3; 
f_amb = 1/(2 .* PRT);


p0 = 0*pi/180; % start angle
p1 = 360*pi/180; % end angle

N_BW = 1;  % Number of beam widths to integrate 

M = round((p1 - p0)/(BW * N_BW)); % Number of azimuth points 



hs = N_BW * round(BW/Omega/PRT); % hits per scan -> Sweeps in one beamwidth

N = hs * M; % Total number of points in time axis

th = linspace(p0, p1, N); % All the angles 
% th = pi/3;
% th = 0;

phi = linspace(th(1), th(end), M); % Angle of the sectors


t1 = 0:PRT:(N - 1)*PRT; % Time axis 

vel_axis = linspace(-f_amb, f_amb, N); % velocity axis for the entire rotation

% mu = linspace(4, 6, 100) .* 2/lambda; % Mean Doppler
[~, idx1] = min(abs(vel_axis - 3 .* 2 / lambda));
[~, idx2] = min(abs(vel_axis - 5 .* 2 / lambda));
mu = vel_axis(idx1:idx2);

% ths = 0:PRT:(hs-1)*PRT;
% 
% t = repmat(ths, [1, M]);
% % 
s_ = zeros(1, length(t1));
s_unmodulated = zeros(1, length(t1));

for i = 1:length(mu)
    ph_ = (2 * pi * mu(i) .* t1);
    s_unmodulated = s_unmodulated + exp(1j .* ph_);
    s_ = s_ + (exp(1j .* ph_ .* (sin(eps + Omega .* t1)./(eps + Omega .* t1))));
end

% s_ = smoothdata(s_);
% s_ = exp(1j .* ph_);
% 

% 
% [~, idx1] = (min(abs(vel_axis - mu)));
% s_f = dirac(vel_axis - vel_axis(idx1)); 
% idx = s_f == Inf;
% s_f(idx) = 1;
% 
% s_ = ifft(ifftshift(sqrt(s_f)));
s_f = fftshift(fft(s_));
figure; plot(vel_axis*lambda/2, db(abs(s_f)), 'LineWidth', 2); grid on; 

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectrum function'); grid on;

s_f_unmodulated = fftshift(fft(hamming(1, N) .* s_unmodulated));
figure; plot(vel_axis*lambda/2, (abs(s_f_unmodulated)), 'LineWidth', 2); grid on; 

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectrum function'); grid on;


s_man = s_;
% s_man = (exp(1j .* unwrap(angle(s_)) .* (eps + sin(th)./(eps + Omega .* t1)))); % manipulated signal with phase correction

% vel_axis = linspace(-v_amb, v_amb, N); % velocity axis for the entire rotation

figure; plot(vel_axis*lambda/2, db(abs(fftshift(fft(s_man))))); % Spectrum of cosine modulated signal
xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spectrum in dB for entire rotation', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectrum function'); grid on;

vel_axis_hs = linspace(-v_amb, v_amb, hs); % velocity axis for one beamwidth

diff_v = diff(vel_axis_hs); 
dv = diff_v(1); % velocity resolution

[S1, F1, Ti1, P1] = spectrogram(s_man, hs, 0, hs, N*Omega/(2*pi));


S1 = S1./max(S1(:));
S1_norm = fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));

F1_doppler = F1 - N/2;

V_doppler = lambda .* F1_doppler ./ 2;

figure;
surface(vel_axis_hs, phi * 180/pi, S1_norm_db); shading flat;
colormap('jet');
colorbar;

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title('Doppler Spectrum');

%% Caluclation of Mean Doppler velocity 


for k = 1:size(S1_norm, 1)
   
    PT = sum(abs(S1_norm(k, :)).^2 .* dv);
    v_mean_i(k, :) = vel_axis_hs .* abs(S1_norm(k, :)).^2 .* dv;
    v_mean(k) = 1./PT .* sum(v_mean_i(k, :));
    
end

txt = ['Angle integration = ', num2str(N_BW*BW_deg), ' [deg]'];
figure(101); hold on; plot(phi * 180/pi, v_mean, 'LineWidth', 2, 'DisplayName', txt); grid on;
legend;
title('Mean Doppler velocity')
xlabel('Angle [deg]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('v_{mean} [m s^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');


