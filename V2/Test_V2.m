clear;
close all;
n = 1024;
N = 128;
PRT = 1e-3;

lambda = 0.03;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
mu = 5;
sigma = 0.2;



phi = pi;

beta_wind = 0;


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
