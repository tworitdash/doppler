%% FFT examples with Gaussian Pulse


clear;
close all;
n = 1024;
N = 1024;
PRT = 1e-3;

lambda = 0.03;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
mu = 5;
sigma = 0.2;



phi = pi;

beta_wind = 0;


[sig, sig_f, sig_f_full, X, T] = DS_simulatorV2(10^(30/10), 1, mu, sigma, n, v_amb, N);
% sig_f = fft(sig(1:256), 256);
% sig_f(end+1) = sig_f(1);
% sig_f2 = fftshift(sig_f);
% 
% sig_fc = fft(conj(sig(1:256)), 256);
% sig_fc(end+1) = sig_fc(1);
% sig_f2c = fftshift(sig_fc);
% 
% % sig = exp(1j .* 2 .* pi .* 2 .* mu / lambda .* (1:N) .* PRT);
% % sig_a = abs(sig) .* exp(1j .* unwrap(angle(sig)) .* cos(beta_wind - phi));
% 
% figure; plot(db(abs(fftshift(fft(sig(1:256))))),'-o');
% figure; plot(db(abs(sig_f2)), '-o'); hold on; plot(flip(db(abs(sig_f2c))));

% hold on; plot(circshift(flip(db(abs(fftshift(fft(conj(sig(1:256))))))), 1));
% hold on; plot(((db(abs(fftshift(fft(conj(sig_a(1:256)))))))));

[sig_with_az, sig_with_az_f] = DS_simulatorV2_with_az(10^(30/10), 1, mu, sigma, n, v_amb, N, phi, beta_wind, X, T);

% sig_ = exp(1j .* 2 .* pi .* 2 .* mu / lambda .* (1:N) .* PRT);
% sig_f = 1./sqrt(N) .* fftshift(fft(sig_, N));
% % 
% sig = ifft(fftshift(sqrt(N) .* sig_f)); % .* exp(1j .* 2 .* pi .* rand(1, N))));
% sig_with_az = exp(1j .* 2 .* pi .* 2 .* (mu .* cos(beta_wind - phi) ) / lambda .* (1:N) .* PRT);

% figure; histogram(real(sig), 1000);
% figure; histogram((sig_f), 1000); 
% figure; histogram((sig_f).^2, 1000);


vel_axis_full = linspace(-v_amb, v_amb, n);
vel_axis = linspace(-v_amb, v_amb, N);


% figure; plot(vel_axis_full, db((abs(sig_f_full).^2))/2); grid on;

% v_content = mean((angle(sig)) .* lambda ./ (4 .* pi .* (1:N) .* PRT)); 
% 
% sig_a = abs(sig) .* exp(1j .* 4 .* pi ./ lambda .* v_content .* cos(beta_wind - phi) .* (1:N) .* PRT);

sig_a = abs(sig) .* exp(1j .* unwrap(angle(sig)) .* cos(beta_wind - phi));

% sig_a = [sig_a1(2:end) sig_a1(1)];

% sig_a = sig .* exp(1j .* 4 .* pi/lambda .* mu .* (1 - cos(beta_wind -  phi)) .* (1:N));

figure; plot(angle(sig_a) .* 180/pi); hold on; plot(angle(sig_with_az) .* 180/pi); legend('my process', 'embedded already')


sig_doppler = 1./sqrt(N) .* abs(fftshift(fft(sig_a, N)));

% sig_with_az_f = 1./sqrt(N) .* abs(fftshift(fft(sig_with_az, N)));

% sig_doppler_max = 1/N .* abs(fftshift(fft(sig, N))).^2;

figure; plot(vel_axis, flip(db(abs(sig_doppler))), '-o'); hold on;  
plot(vel_axis_full, (db(abs(sig_f_full))));hold on;  
plot(vel_axis, flip(db(abs(sig_with_az_f))));  grid on; legend('with cosine', 'actual', 'with cosine embedded')



sig_doppler_req = 1./sqrt(N) .* abs ((fft(sig_a, N)));
% figure; plot(vel_axis, db(abs(sig_doppler_req)), '-o'); title('without fftshift')


if cos(beta_wind - phi) < 0
    vel_axis_proc = (vel_axis(1:N/2));
     Signal_dop_required = sig_doppler_req(N/2+1:N);
%                     vel_axis_proc(m, i).axis = vel_axis(m, i).axis(hits_scan_(m, i)/2:hits_scan_(m, i));
else
    vel_axis_proc = vel_axis([N/2+1:N]);
    Signal_dop_required = sig_doppler_req(1:N/2);
end


PRT = 1e-3; lambda = 3e-2;

del_v = lambda/(2*N*PRT);
                
figure; plot(vel_axis_proc, db(abs(squeeze(Signal_dop_required))), 'o', 'color', 'k');
                
PT_integrand = abs(Signal_dop_required).^2 .* del_v;
PT = sum(PT_integrand); % Total power of the Doppler Spectrum

v_mean_integrand = vel_axis_proc .* abs(Signal_dop_required).^2 .* del_v;
v_mean_l = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

v_spread_integrand = (vel_axis_proc - v_mean_l).^2 .* abs(Signal_dop_required).^2 .* del_v;
v_spread_l = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
                



PT_integrand = abs(sig_doppler).^2 .* del_v;
PT = sum(PT_integrand); % Total power of the Doppler Spectrum

v_mean_integrand = vel_axis .* abs(sig_doppler).^2 .* del_v;
v_mean = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

v_spread_integrand = (vel_axis - v_mean_l).^2 .* abs(sig_doppler).^2 .* del_v;
v_spread = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width

%%
phi = pi;
beta_wind = 0;

sig = exp(1j .* 2 .* pi .* 2 .* mu / lambda .* (1:N) .* PRT);
sig_a = abs(sig) .* exp(1j .* unwrap(angle(sig)) .* cos(beta_wind - phi));
sig_c = exp(1j .* 2 .* pi .* 2 .* (mu .* cos(beta_wind - phi)) ./ lambda .* (1:N) .* PRT);
N = 256;
sig_f = 1./sqrt(N) .* fftshift(fft(sig, N));
sig_af = 1./sqrt(N) .* fftshift(fft(sig_a, N));
sig_cf = 1./sqrt(N) .* fftshift(fft(sig_c, N));

figure; plot(db(abs(sig_f))); hold on; plot(flip(db(abs(sig_af)))); hold on; plot(flip(db(abs(sig_cf))))

