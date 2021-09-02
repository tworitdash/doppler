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
mu = 7.5/3 .* 2./lambda;

p = 1; % skew of the spectrum
wts = linspace(0.03, 0.28, 10);
% wts = 0.28;
sigma = wts/(p .* PRT);



v_amb = lambda/(4 * PRT) .* 2./lambda;

SNR_db = 60;
SNR = 10^(SNR_db/10);
%% Gathering time domain data from a Gaussian spectrum

for m = 1:length(sigma)

[data, data_f, data_f_Sig, X, Theta, P_wNoise, P] = DS_simulatorV3(SNR, m0, mu, sigma(m), N, v_amb, p);
data = data.';

if mod(N, 2) == 0
    vel_axis = linspace(-N/2, N/2-1, N)./N .* 2 .* v_amb;
else
    vel_axis = linspace(-v_amb, v_amb, N);
end


figure; plot(vel_axis/v_amb, db(P_wNoise./max(P_wNoise))./2, 'LineWidth', 2); grid on;
figure; plot(vel_axis/v_amb, (-P./max(-P)), 'LineWidth', 2); grid on;




%% Pulse pair and Poly Pulse pair

k = 1:1:5;
% k = 1;
% R = zeros(sec, length(k));

for l = 1:length(k)
    R(l) = 1./N .* [data(1:end-k(l)+1)'] * [data(k(l):end)];
    fn(m, l) = 1./(2 .* pi .* PRT .* k(l)) .* angle(R(l));
end

% figure(101); plot(k, fn .* lambda/2, 'LineWidth', 2); grid on; hold on; plot(k, mu .* lambda/2 .* ones(1, length(k)), 'LineWidth', 2);

f_vn(m) = (1/(2 * pi * PRT)) .* angle(sum(abs(R) .* exp(1j .* 2 * pi .* fn(m, :) .* k .* PRT))) .* lambda/2;
    
f_vun(m) = (1/(2 * pi * PRT)) .* angle(sum(abs(R).^(1./k) .* exp(1j .* 2 * pi .* fn(m, :) .* k .* PRT))) .* lambda/2;

f_av(m) = 1/length(k) .* sum(fn(m, :)) .* lambda/2;

end

figure; plot(wts,(fn(:, 2)));
hold on; plot(wts, (fn(:, 4)));
hold on;  plot(wts, (fn(:, 5)));

grid on;

figure; plot(wts, f_vn); hold on;
plot(wts, f_vun); hold on;
plot(wts, f_av); hold on; grid on;

% hold on; plot(k, f_vun .* lambda/2, 'o'); 
% hold on; plot(k, f_av .* lambda/2, 'LineWidth', 2);
% hold on; plot(k, mu .* lambda/2 .* ones(1, length(k)), '*');
