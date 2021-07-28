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
% mu_r = 1./2;
mu_r = linspace(-1, 1, 10);
mu = v_amb .* mu_r;

p = 1.5; % skew of the spectrum
% wts = linspace(0, 0.14, 10);
wts = 0.0625;
% sigma = wts/(p .* PRT);
sigma = wts / (2 .* PRT);

% SNR_db = linspace(0, 40, 4);
SNR_db = 40;
SNR = 10.^(SNR_db/10);
%% Gathering time domain data from a Gaussian spectrum
for s = 1:length(SNR)
for m = 1:length(mu)

[data, data_f, data_f_Sig, X, Theta, P_wNoise, P] = DS_simulatorV3(SNR(s), m0, mu(m), sigma, N, v_amb, p);
data = data.';

if mod(N, 2) == 0
    vel_axis = linspace(-N/2, N/2-1, N)./N .* 2 .* v_amb;
else
    vel_axis = linspace(-v_amb, v_amb, N);
end

dv = vel_axis(2) - vel_axis(1);
% 
% figure; plot(vel_axis/v_amb, db(P_wNoise./max(P_wNoise))./2, 'LineWidth', 2); grid on;
% figure; plot(vel_axis/v_amb, (-P./max(-P)), 'LineWidth', 2); grid on;


I = real(data);
Q = imag(data);

%% DFT 

% S_FFT = 1/sqrt(N) .* fftshift(fft(data));
% 
% PT = sum(abs(S_FFT).^2 .* dv);
% 
% M1_FFT(m) = sum(vel_axis.' .* abs(S_FFT).^2 .* dv)./PT;
% 
% M2_FFT(m) = sqrt(sum((vel_axis.' - M1_FFT(m)).^2 .* abs(S_FFT).^2 .* dv) ./ PT);


%% Pulse pair and Poly Pulse pair

R = 1./N .* data(2:end).' * (data(1:end-1)').';

R0 = 1./N .* data.' * (data').';

M1_PP(m) = 1./(2 .* pi .* PRT) .* angle(R);

M2_PP(m) = 1/(pi .* PRT .* sqrt(2)) .* sqrt(1 - abs(R)./R0);


%% Vector Phase Change
num = 0;
dum = 0;
for i = 1:N - 1
    num = num + sin(angle(data(i + 1)) - angle(data(i)));
    dum = dum + cos(angle(data(i + 1)) - angle(data(i)));
end

M1_VPC(m) = (1./(2 .* pi .* PRT)) .* angle(dum + 1j .* num);


%% Poly Pulse Pair

k = 1:1:round(N/2)- 1;

for l = 1:length(k)
    R_pppn(l) = 1./N .* data(k(l)+1:end).' * (data(1:end-k(l))').';
    f_pppn(l) = 1./k(l) .* 1./(2 .* pi .* PRT) .* angle(R_pppn(l));
end


%% VPP-N

f_ppp_vppn(m) = 1./(2 .* pi .* PRT) .* angle(sum(R_pppn));

end



% 
% txt = ['Input SNR = ', num2str(SNR_db(s)), ' dB'];
% figure(1); hold on; plot(wts, 2 .* PRT .* abs(M2_PP), 'DisplayName', txt, 'LineWidth', 2); grid on; 
txt = ['Input SNR = ', num2str(SNR_db(s)), ' dB']; figure(101); hold on;
plot(mu_r, f_ppp_vpn./v_amb, 'DisplayName', txt, 'LineWidth', 2); grid on;
end

% xlabel('2 \sigma_{fTrue} PRT', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('2 \sigma_{fPP} PRT', 'FontSize', 12, 'FontWeight', 'bold');
% title('Standard Deviation PP Estimate');
% legend; 


xlabel('\mu_{fTrue}/{f_{amb}}', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\mu_{fPPPVPPN}/{f_{amb}}', 'FontSize', 12, 'FontWeight', 'bold');
title('Mean Doppler PPP Estimate with VPP-N');
legend; 



