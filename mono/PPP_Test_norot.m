clear;
close all;
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
% phi_axis = linspace(phi_0, Omega*t(end), sec);
% phi_axis = zeros(size(phi_axis));
m0 = 1;
% mu = 2 .* 2./lambda;

% mu = 0;

p = 2; 
g = 4 * sqrt(2)/sqrt(pi) * [(p^(-2) + 1)^(-1.5) - (p^2 + 1)^(-1.5)]; % Skewness

wts = linspace(0.03, 0.28, 10);

sigma = wts/(PRT);
v_amb = lambda/(4 * PRT) .* 2./lambda;


% mu = linspace(-v_amb, v_amb, 10);
mu = 1 * 2/lambda;
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

% % 
% figure; plot(vel_axis/v_amb, db(P_wNoise./max(P_wNoise))./2, 'LineWidth', 2); grid on;
% figure; plot(vel_axis/v_amb, (-P./max(-P)), 'LineWidth', 2); grid on;
% 
% figure; plot(vel_axis*2/lambda, db(abs(fftshift(fft(data)))));


%% Pulse pair and Poly Pulse pair

k = 1:1:10;


for l = 1:length(k)
    lag = k(l);
    
    Rel(m, l) = 1/(N) * (data(1:end-lag))' * data(lag+1:end);
    fn(m, l) = 1/(2 * pi * lag * PRT) * angle(Rel(m, l));
    fn_err(m, l) = abs(fn(m, l) - mu)./(mu) * 100;
end


f_vppn(m) = 1/(2 * pi * PRT) * angle(sum(abs(Rel(m, :)) .* exp(1j * 2 * pi * fn(m, :) * PRT)));



end



figure; plot(sigma*lambda/2,(fn_err(:, 1)));
hold on; plot(sigma*lambda/2,(fn_err(:, 3)));
hold on;  plot(sigma*lambda/2, (fn_err(:, 5)));

grid on;



figure; plot(sigma*lambda/2,lambda/2*(fn(:, 1)));
hold on; plot(sigma*lambda/2,lambda/2*(fn(:, 3)));
hold on;  plot(sigma*lambda/2, lambda/2*(fn(:, 5)));

grid on;
% hold on; plot(sigma * lambda/2, f_vppn * lambda/2); hold on;


% plot(wts, f_vun); hold on;
% plot(wts, f_av); hold on; grid on;

% hold on; plot(k, f_vun .* lambda/2, 'o'); 
% hold on; plot(k, f_av .* lambda/2, 'LineWidth', 2);
% hold on; plot(k, mu .* lambda/2 .* ones(1, length(k)), '*');

