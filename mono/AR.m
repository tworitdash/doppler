%% AR modeling for Doppler estimation 

clear;
close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = 6;

lambda = 0.03;

N = 2048;

t = 0:PRT:(N - 1)*PRT;

m0 = 1;

p = 1; 
g = 4 * sqrt(2)/sqrt(pi) * [(p^(-2) + 1)^(-1.5) - (p^2 + 1)^(-1.5)]; % Skewness

% wts = linspace(0.03, 0.28, 10);
wts = 0.1;

sigma = wts/(PRT);
v_amb = lambda/(4 * PRT) .* 2./lambda;

mu = 1 * 2/lambda;
SNR_db = 40;
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

end
eta = sqrt(0.1) .* rand(1, N);
I_ = real(data) + eta;
Q_ = imag(data) + eta;
data = I_ + 1j .* Q_;
a = zeros(1,  N);
Z = zeros(1, N);
gamma = zeros(1, N);


a(1) = 0;
Z(1) = 0;
gamma(1) = abs(data(1))^(-2);
P(1) = 1;

for i = 1:N-1
    Q(i) = sum(abs(data(1:i)).^2)./(sum(eta(1:i).^2 - 1));
    P(i + 1) = ((Q(i) .* (1 - abs(a(i)).^2)) + P(i) .* abs(a(i)).^2)...
        ./(1 + ((Q(i) .* (1 - abs(a(i)).^2)) + P(i) .* abs(a(i)).^2));
    Z(i + 1) = a(i) * Z(i) + P(i + 1) .* (data(i + 1) - a(i) .* Z(i));
    gamma(i + 1) = gamma(i)./(1 + gamma(i) * abs(Z(i)).^2);
    a(i+1) = a(i) + gamma(i + 1) .* (Z(i))' * (data(i + 1) - a(i) * Z(i));
end

figure; plot(angle(a));
