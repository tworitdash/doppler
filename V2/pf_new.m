clear;
% close all;


lambda = 0.03;

PRT = 1e-3;
mu = 5; % Mean Doppler
sigma = 0.2;
v_amb = 7.5;

N = 60000; % Total number of points in time axis

t1 = 0:PRT:(N - 1)*PRT; % Time axis
vel_axis = linspace(-v_amb, v_amb, N); % velocity axis for the entire rotation

% vel_axis_hs = linspace(-v_amb, v_amb, hs); % velocity axis for one beamwidth
% t1 = -N/2*PRT:PRT:(N/2 - 1)*PRT; % Time axis
% s_analyt_o = 1./sqrt(2*pi*(sigma).^2) .* exp(-(vel_axis - mu).^2./(2*(sigma).^2));

% SNR = 10^(30/10);
% X = rand(1, N);
% Theta = 2 .* pi * rand(1, N);


S_ = gaussmf(vel_axis, [sigma, mu]);

% Noise = sum(S_) ./ (N .* SNR);
% s_analyt_o = -(S_ + Noise) .* log(X);

s_analyt_o = S_;
s_num = ifft(ifftshift(sqrt(s_analyt_o) .* exp(1j .* 2 * pi * rand(1, N)))); % Numerical IFFT


s_analyt = ((2/pi)^(1/4) * sqrt(sigma * 4 * pi / lambda) .* exp(-t1 .* ((sigma * 4 * pi / lambda)^2 .* t1 - 1j .* mu .* 4 * pi / lambda))); % analytical IFFT



figure;
plot(t1, abs((s_analyt))/max(abs(s_analyt)), 'LineWidth', 2);

hold on; plot(t1, abs((s_num))/max(abs(s_num)), '*');
grid on;

xlim([-10 80]);

legend({'Wolfram Alpha Expression', 'MATLAB numerical'})


figure; 

plot(vel_axis, abs(fftshift(fft(s_analyt)))/max(abs(fft(s_analyt)))); hold on;
plot(vel_axis, abs(fftshift(fft(s_num)))/max(abs(fft(s_num))));


legend({'Wolfram Alpha Expression', 'MATLAB numerical'})