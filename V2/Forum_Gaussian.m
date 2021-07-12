clear; 
% close all;


lambda = 0.03;


mu = 4 * 2 / lambda; % Mean Doppler
sigma = 0.1 * 2 / lambda;


PRT = 1e-3; 
f_amb = 1/(2 * PRT);

hs = 300;
N = hs * 200; % Total number of points in time axis


f = linspace(-f_amb, f_amb, N); % velocity axis for the entire rotation



S_ = gaussmf(f, [sigma, mu]);


s_num = ifft(ifftshift(sqrt(S_)));

t1 = 0:PRT:(N - 1) * PRT;

th = linspace(0, 2*pi, N);

Omega = (1.8 * pi/180)/(hs .* PRT);

s_analyt = (2/pi)^(1/4) * sqrt(sigma.*cos(th)) .* exp(-(sigma.*cos(th)).^2.*t1.^2) .* exp(1j .* 2 .* pi .* mu .* t1 .* (eps + sin(Omega .* t1))./(eps + Omega .* t1)); % analytical IFFT 
% s_analyt = hamming(1, N) .* ifftshift((2/pi)^(1/4) * sqrt(sigma) .* exp(sigma^2*t1.^2) .* exp(1j .* 2 .* pi .* mu .* t1)); % analytical IFFT s_analyt = 0 .* exp(1j .* 2 .* pi .* mu .* t1);



figure; plot(unwrap(angle(s_num))/(2*pi*mu));
figure; plot(unwrap(angle(s_analyt))/(2*pi*mu));

[S1, F1, Ti1, P1] = spectrogram(s_analyt, hs, 0, hs , N*1/60);


S1 = S1./max(S1(:));
S1_norm = fftshift(S1',2);
S1_norm_db = 20*log(abs(S1_norm));

% F1_doppler = F1 - N/2;
% 
% V_doppler = lambda .* F1_doppler ./ 2;


f_axis_hs = linspace(-f_amb*lambda/2, f_amb*lambda/2, hs);


figure;
surface(f_axis_hs, Ti1, S1_norm_db); shading flat;
colormap('jet');
colorbar;

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
title('Doppler Spectrum');


%% 




[~, idx1] = (min(abs(f - mu)));
S_Dirac_f = dirac(f - f(idx1)); 
idx = S_Dirac_f == Inf;
S_Dirac_f(idx) = 1;

s_num_dirac = ifft(ifftshift(sqrt(S_Dirac_f)));

figure; plot(unwrap(angle(s_num_dirac))/(2*pi*mu));