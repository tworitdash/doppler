clear; 
close all;


lambda = 0.03;


mu = 4; % Mean Doppler
Omega_rpm = 1; % in RPM
Omega = 2*pi/60 * Omega_rpm; % In rad/s
BW_deg = 1.8; % beam width in degree
BW = BW_deg * pi/180; % beam width in radian

v_amb = 7.5;

PRT = 1e-3; 
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


% t1_ = 0:PRT:(N - 1)*PRT; % Time axis 
% t1 = 0:N-1;


% ths = 0:PRT:(hs-1)*PRT;
% 
% t = repmat(ths, [1, M]);
% % 

% ph_ = (4 * pi * mu / lambda .* t1);
% % s_ = (exp(1j .* ph_ .* (sin(eps + Omega .* t1)./(eps + Omega .* t1) - 1)));
% s_ = exp(1j .* ph_);
% 
vel_axis = linspace(-v_amb, v_amb, N); % velocity axis for the entire rotation
% 
[~, idx1] = (min(abs(vel_axis - mu)));
s_f = dirac(vel_axis - vel_axis(idx1)); 
idx = s_f == Inf;
s_f(idx) = 1;

s_ = ifft(ifftshift(sqrt(s_f)));

s_man = (exp(1j .* unwrap(angle(s_)) .* (eps + sin(th)./(eps + Omega .* t1_)))); % manipulated signal with phase correction

s2 = ifftshift(exp(1j .* 4 * pi / lambda * mu * (eps + sin(th)./(eps + Omega))));

figure; plot(t1, real(s_man)); hold on; plot(t1, real(s2));

figure; plot(vel_axis, db(abs(fftshift(fft(s_man))))); hold on; plot(vel_axis, db(abs(fftshift(fft(s2)))), '*');



