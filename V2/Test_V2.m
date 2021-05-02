
%% Radar specification
clear;
close all;
n = 1024;
N = 8;
PRT = 1e-3;

lambda = 0.03;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
mu = 5;
sigma = 0.2;



phi = pi/4;

beta_wind = 0;

%% Signal and spectrum generation

[sig, sig_f, sig_f_full, sig_full, X, T] = DS_simulatorV2(10^(30/10), 1, mu, sigma, n, v_amb, N);

[sig_with_az, sig_with_az_f] = DS_simulatorV2_with_az(10^(30/10), 1, mu, sigma, n, v_amb, N, phi, beta_wind, X, T);


vel_axis_full = linspace(-v_amb, v_amb, n);
vel_axis = linspace(-v_amb, v_amb, N);

%% First procedure  (Time domain sampling)

idx = 1:n;
idxq = linspace(min(idx), max(idx), N);

sig_resampled = interp1(idx, sig_full, idxq, 'nearest');
time_axis_full = eps:PRT:(n-1)*PRT;
time_axis_resampled = linspace(min(time_axis_full), max(time_axis_full), N);
 
figure; plot(time_axis_full, real(sig_full), 'LineWidth', 2); hold on; plot(time_axis_resampled, real(sig_resampled), '*');legend({'Original Signal real part', 'Signal real part after resampling'});
grid on;
sig_f_resampled = 1/sqrt(N) .* fftshift(fft(sig_resampled));

figure; plot(vel_axis_full, db(abs(sig_f_full)), 'LineWidth', 2); hold on; plot(vel_axis, db(abs(sig_f_resampled)), 'LineWidth', 2); legend({'Original Spectrum', 'Spectrum after resampling'});
grid on;

%% Second procedure (Frequency domain sampling) indpendent of direction information

figure; plot(vel_axis, db(abs(sig_f)), 'LineWidth', 2); title('Frequency domain resampling'); grid on;

figure; plot(time_axis_full, real(sig_full), 'LineWidth', 2); hold on; plot(time_axis_resampled, real(sig), '*');legend({'Original Signal real part', 'Signal real part after resampling'});
grid on;
sig_doppler = 1/sqrt(N) .* fftshift(fft(sig));

figure; plot(vel_axis_full, db(abs(sig_f_full)), 'LineWidth', 2); hold on; plot(vel_axis, db(abs(sig_doppler)), 'LineWidth', 2); legend({'Original Spectrum', 'Spectrum after resampling'});
grid on;

%% Symmetry Problem of frequency domain sampling independent of direction information 
Nfft = N + 1;
vel_axis_fft = linspace(-v_amb, v_amb, Nfft);
vel_axis = linspace(-v_amb, v_amb, N);

sig_a = abs(sig) .* exp(1j .* unwrap(angle(sig)) .* cos(beta_wind - phi));

sig_a_doppler = 1/sqrt(Nfft) .* fftshift(fft(sig_a, Nfft)); 

figure; plot(vel_axis, db(abs(sig_doppler)), 'LineWidth', 2); hold on; plot(vel_axis, flip(db(abs(sig_with_az_f))), 'LineWidth', 2) 
legend({'Spectrum after resampling in frequency domain', 'Spectrum after using \Phi = \pi and flipped'}, 'FontSize', 16);
grid on;

%% Monochromatic signal 
% v_max = 7.5;
% v = 5;
% N = 128;
% PRT = 1e-3;
% vel_axis = linspace(-v_max, v_max, N);
% x = (exp(1j .* 4 .* pi .* v ./ lambda .* (1:N) .* PRT));
% x_conj = (conj(exp(1j .* 4 .* pi .* v ./ lambda .* (1:N) .* PRT)));
% x_fft = 1/sqrt(N) .* fftshift(fft(x, N));
% x_conj_fft = 1/sqrt(N) .* fftshift(fft(x_conj, N));
% figure; plot(vel_axis, db(abs(x_fft)), 'LineWidth', 2); hold on; plot(vel_axis, (db(abs(x_conj_fft))), 'LineWidth', 2);grid on;
% 
% legend({'Spectrum', 'Conjugate Spectrum'}, 'FontSize', 16);

%% 
v_max = 7.5;
v = 5;
N = 8;
Nfft = 128;
PRT = 1e-3;
lambda=0.03;
Axis=((0:Nfft-1) -ceil((Nfft-1)/2))/Nfft;
t_axis=1:N;
vel_axis_ =2*v_max*Axis; %<--- Changed
% x = ifftshift( exp(1j .* 4 .* pi .* v ./ lambda .* t_axis .* PRT) );
% % x_conj =  conj(x) ;
% x_conj = ifftshift( exp(1j .* 4 .* pi .* v .* cos(beta_wind - phi) ./ lambda .* t_axis .* PRT) );
sig_a1 = ifftshift(sig_a);
sig_a_doppler = 1/sqrt(Nfft) .* (fft(sig_a1, Nfft)); 
% flipit=@(z) z([1,end:-1:2]) ; %<--- Added
% x_fft = 1/sqrt(N) .* fft(x, N);
% x_conj_fft = 1/sqrt(N) .* fft(x_conj, N); %<--- Changed
% [x_max,idx]=max(db(abs(x_fft)));
%   x_max
% x_max = 19.4222
%   vel_max=vel_axis(idx)
% vel_max = -2.4609
%  [xconj_max,idx]=max(db(abs(x_conj_fft)));
%   xconj_max
% xconj_max = 19.4222
%   velconj_max=vel_axis(idx) 
% velconj_max = 2.4609
  
  
% figure; plot(vel_axis, fftshift(db(abs(x_fft))), 'LineWidth', 2); 
% hold on; 
% plot(vel_axis, fftshift(db(abs(x_conj_fft))), 'LineWidth', 2);
% hold on;
% plot(vel_axis, (db(abs(sig_with_az_f))), 'LineWidth', 2) 
% hold off
% grid on;

figure; plot(vel_axis_full, db(abs(sig_f_full)), 'LineWidth', 2);
hold on;
plot(vel_axis_, fftshift(db(abs(sig_a_doppler))), 'LineWidth', 2);
hold on;
plot(vel_axis, (db(abs(sig_with_az_f))), 'LineWidth', 2) 
hold off
grid on;

%%

v_max = 7.5;
v = 5;
N = 128;
Nfft = 128;
PRT = 1e-3;
lambda=0.03;
Axis=((0:Nfft-1) -ceil((Nfft-1)/2))/Nfft;
t_axis=1:N;
vel_axis =2*v_max*Axis; %<--- Changed
x = ifftshift( exp(1j .* 4 .* pi .* v ./ lambda .* t_axis .* PRT) );
x_conj =  conj(x) ;
% x_conj = ifftshift( exp(1j .* 4 .* pi .* v .* cos(beta_wind - phi) ./ lambda .* t_axis .* PRT) );
% sig_a1 = ifftshift(sig_a);
% sig_a_doppler = 1/sqrt(N) .* (fft(sig_a1, N)); 
flipit=@(z) z([1,end:-1:2]) ; %<--- Added
x_fft = 1/sqrt(Nfft) .* fft(x, Nfft);
x_conj_fft = 1/sqrt(Nfft) .* fft(x_conj, Nfft); %<--- Changed
[x_max,idx]=max(db(abs(x_fft)));
%   x_max
% x_max = 19.4222
  vel_max=vel_axis(idx)
% vel_max = -2.4609
 [xconj_max,idx]=max(db(abs(x_conj_fft)));
%   xconj_max
% xconj_max = 19.4222
  velconj_max=vel_axis(idx) 
% velconj_max = 2.4609
  
  
figure; plot(vel_axis, fftshift(db(abs(x_fft))), 'LineWidth', 2); 
hold on; 
plot(vel_axis, fftshift(db(abs(x_conj_fft))), 'LineWidth', 2);
hold on;
% plot(vel_axis, flipit(fftshift(db(abs(x_conj_fft)))), '*', 'LineWidth', 2);
% plot(vel_axis, (db(abs(sig_with_az_f))), 'LineWidth', 2) 
% hold off
grid on;

% figure; plot(vel_axis_full, db(abs(sig_f_full)), 'LineWidth', 2);
% hold on;
% plot(vel_axis, fftshift(db(abs(sig_a_doppler))), 'LineWidth', 2);
% hold on;
% plot(vel_axis, (db(abs(sig_with_az_f))), 'LineWidth', 2) 
% hold off
% grid on;