clear;
% close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = 6;

lambda = 0.03;

BW = BW_deg * pi/180;
Omega = Omega_rpm * 2*pi/60; % rotation speed in rad/s
Td = BW/Omega;
hs = round(Td/PRT);
phi_0 = phi_0_deg * pi/180;

sec = round((n_rot*2*pi)/BW);
% N = sec * hs;
N = 2048;

t = 0:PRT:(N - 1)*PRT;
% phi_axis = phi_0:BW:phi_0+Omega*t(end);
phi_axis = linspace(phi_0, Omega*t(end), sec);
% phi_axis = zeros(size(phi_axis));
m0 = 1;
v_amb = lambda/(4 * PRT) .* 2./lambda;
mu_r = linspace(-1, 1, 100);
mu = v_amb .* mu_r;

p = 1; % skew of the spectrum
% wts = linspace(0.03, 0.28, 10);
% wts = 0.28;
% sigma = wts/(p .* PRT);
sigma = 0.0625 / (2 .* PRT);

SNR_db = 30;
SNR = 10^(SNR_db/10);
%% Gathering time domain data from a Gaussian spectrum

for m = 1:length(mu)

[data, data_f, data_f_Sig, X, Theta, P_wNoise, P] = DS_simulatorV3(SNR, m0, mu(m), sigma, N, v_amb, p);
data = data.';

if mod(N, 2) == 0
    vel_axis = linspace(-N/2, N/2-1, N)./N .* 2 .* v_amb;
else
    vel_axis = linspace(-v_amb, v_amb, N);
end


% figure; plot(vel_axis/v_amb, db(P_wNoise./max(P_wNoise))./2, 'LineWidth', 2); grid on;
% figure; plot(vel_axis/v_amb, (-P./max(-P)), 'LineWidth', 2); grid on;


I = real(data);
Q = imag(data);
% figure(102); hold on; plot(vel_axis/v_amb, db(1./N .* abs(fftshift(fft(data)))), 'LineWidth', 2); grid on;
% 
% figure(m); plot(t, I, 'LineWidth', 2); grid on;
% hold on; plot(t, Q, 'LineWidth', 2); grid on;

Num = 0;
Dum = 0;
for n = 2:N
    Num = Num + I(n - 1) .* (Q(n) - Q(n - 1)) + Q(n - 1) .* (I(n) - I(n - 1));
    Dum = Dum + I(n - 1).^2 + Q(n - 1).^2; 
end

%% Pulse pair and Poly Pulse pair

fTDC(m) = Num./Dum;
end

figure; plot(mu_r, fTDC/v_amb);




