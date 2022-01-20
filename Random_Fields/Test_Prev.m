clear;
close all;

m_ = load('../mono/markers.mat'); 
markers = m_.markers;

%% Simulation Inputs

PRT = 1e-3;
lambda = 0.03;

v_amb = lambda/(4 * PRT);

mu = 0;
sigma = eps;
N_rot = 1;
Sweep = 1000;

n = N_rot * Sweep;

SNR_db = 30;
SNR = 10^(SNR_db/10);

m0 = 1;

%% Simulator Output
[df, d, ~, ~, vel_axis] = DS(SNR, m0, mu, sigma, n, v_amb);

% txt = ['Original Spectrum'];
% 
% dtext = ['SNR = ', num2str(SNR_db), ' dB', ', \mu = ', num2str(mu), ' m/s', ' \sigma = ', num2str(sigma), ' m/s'];
% xl = 'velocity [m/s]';
% 
% f = figure(105); hold on; f.Position = [10 10 1000 1000];
% color = 'k';
% yl =  ['Spectrum [dB]'];
% plott(vel_axis, db(abs(df)), xl, yl, txt, 2, dtext, color, markers(1));

%% Process

Nfft = 5;

df_process = 1/sqrt(Nfft) .* fftshift(fft(d(1:Nfft), Nfft));
% df_process = 1/sqrt(n) .* (fft(d, Nfft));
if mod(Nfft, 2) == 0
    
    axis = linspace(-Nfft/2, Nfft/2-1, Nfft)/(Nfft);
    vel_axis_n = 2 * v_amb * axis;
    f_axis_n = linspace(0, 2*v_amb*2/lambda, Nfft);
    
else
    
    vel_axis_n = linspace(-v_amb, v_amb, Nfft);
    f_axis_n = linspace(0, 2*v_amb*2/lambda, Nfft);
end

txt = ['Spectrum'];

dtext = ['SNR = ', num2str(SNR_db), ' dB', ', \mu = ', num2str(mu), ' m/s', ' \sigma = ', num2str(sigma), ' m/s'];
xl = 'frequency [Hz]';

f = figure(105); hold on; f.Position = [10 10 1000 1000];
hold on;
color = 'r';
yl =  ['Processed Spectrum [dB] from simulated data'];
plott2(vel_axis_n, db(abs(df_process)), xl, yl, txt, 2, dtext, color, markers(2))

dv = vel_axis_n(2) - vel_axis_n(1);
PT = sum(abs(df_process).^2  .* dv);
mu = 1./PT .* sum(vel_axis_n .* abs(df_process).^2 .* dv);
sigma = sqrt(1./PT .* sum((vel_axis_n - mu).^2 .* abs(df_process).^2 .* dv));
