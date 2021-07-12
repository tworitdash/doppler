clear; 
close all;


lambda = 0.03;


mu = 4 * 2 / lambda; % Mean Doppler
sigma = 0.1 * 2 / lambda;
Omega_rpm = 1; % in RPM
Omega = 2*pi/60 * Omega_rpm; % In rad/s
BW_deg = 1.8; % beam width in degree
BW = BW_deg * pi/180; % beam width in radian

v_amb = 7.5;

PRT = 1e-3; 
f_amb = 1/(2 * PRT);

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
% t1 = (-N/2)*PRT:PRT:(N/2-1)*PRT;
% t1 = 0:N-1;
vel_axis = linspace(-f_amb, f_amb, N); % velocity axis for the entire rotation

vel_axis_hs = linspace(-v_amb, v_amb, hs); % velocity axis for one beamwidth
% t1 = -N/2*PRT:PRT:(N/2 - 1)*PRT; % Time axis
% s_analyt_o = 1./sqrt(2*pi*(sigma).^2) .* exp(-(vel_axis - mu).^2./(2*(sigma).^2));

SNR = 10^(30/10);
X = rand(1, N);
Theta = 2 .* pi * rand(1, N);


S_ = gaussmf(vel_axis, [sigma, mu]);

Noise = sum(S_) ./ (N .* SNR);

s_analyt_o = -(S_ + Noise) .* log(X);
% s_analyt_o = -S_;
% figure; plot(vel_axis*lambda/2, db(abs(s_analyt_o))/2); grid on; 
% xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Spectrum without modulation/rotation in dB', 'FontSize', 12, 'FontWeight', 'bold');
% title('Doppler Spectrum');
% 

s_analyt = exp(- (2 * pi * sigma .* cos(th)).^2) .* exp(t1 .* (1j .* mu .* 2.*pi .* (eps + sin(th))./(eps + Omega .* t1)));


s_num = ifft(ifftshift(sqrt(s_analyt_o) .* exp(1j .* Theta)));
% 
figure; plot(th*180/pi, abs((s_analyt)));
figure; plot(vel_axis*lambda/2, db(abs(fftshift(fft(s_analyt)))));
% 
% figure; plot(t1, abs(ifftshift(s_num)));

[S1, F1, Ti1, P1] = spectrogram(s_analyt, hs, 0,hs , N*Omega/(2*pi));


S1 = S1./max(S1(:));
S1_norm = fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));

% F1_doppler = F1 - N/2;
% 
% V_doppler = lambda .* F1_doppler ./ 2;

figure;
surface(vel_axis_hs, phi * 180/pi, S1_norm_db); shading flat;
colormap('jet');
colorbar;

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title('Doppler Spectrum');


s_analyt_f = fftshift(fft(s_analyt)); 

figure; plot(vel_axis*lambda/2, db(abs(s_analyt_f)));

% s_man = (2/pi)^(1/4) * sqrt(sigma * 4 * pi /lambda ) .* exp(-t1 .* ((sigma * 4 * pi /lambda)^2  .* t1 - 1j .* mu .* 4 * pi /lambda .* ((eps + sin(th))./(eps + th))));
% s_man = (2/pi)^(1/4) * sqrt(sigma * 4 * pi /lambda ) .* exp(-t1 .* (- 1j .* mu .* 4 * pi /lambda .* ((eps + sin(th))./(eps + th))));
% s_man = s_analyt;
s_man = abs((s_num)) .* exp(1j .* unwrap(angle((s_num))) .* (eps + sin(Omega .* t1))./((eps + Omega .* t1)));
 
% vel_axis = linspace(-v_amb, v_amb, N); % velocity axis for the entire rotation

figure; plot(vel_axis, db(abs(fftshift(fft(s_man))))); % Spectrum of cosine modulated signal
xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spectrum in dB for entire rotation in dB', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectrum function'); grid on;



diff_v = diff(vel_axis_hs); 
dv = diff_v(1); % velocity resolution

[S1, F1, Ti1, P1] = spectrogram(s_man, hs, 0, hs, N*Omega/(2*pi));


S1 = S1./max(S1(:));
S1_norm = fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));

F1_doppler = F1 - N/2;

V_doppler = lambda .* F1_doppler ./ 2;

figure;
imagesc(vel_axis_hs, phi * 180/pi, S1_norm_db); shading flat;
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


