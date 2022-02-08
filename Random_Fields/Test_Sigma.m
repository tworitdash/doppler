%%
clear;
close all;
markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

N = 128;

r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 0.01*pi/180;

r = r0 - dr/2 + dr .* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);

% x0 = 10 + 5 .* rand(1, N);
% % y0 = 10 + 5 .* rand(1, N); 
% y0 = zeros(1, N);

% x0 = r .* cos(ph); y0 = r .* sin(ph);


x0 = linspace(50, 100, N);
y0 = zeros(1, N);
x(1, :) = x0; y(1, :) = y0;
 
%% Plot positions of scatterers
txt = ['Resolution cell'];
dtext = ['Scaterer position at time t = 0'];
xl = 'x[m]';

f = figure(107); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['y[m]'];

marker = markers(2);

plott2(x0, y0, xl, yl, txt, 2, dtext, color, marker)
%%
% figure(101); plot(x0, y0, '*');

A0 = ones(1, N);
D0 = 0;

lambda = 0.03; dt = 1e-3;

z(1) = sum(A0 .* exp(1j .* 4 * pi / lambda * x0));

SNR_db = 60; 
SNR = 10^(SNR_db/10);

Vspread = linspace(eps, 3, 100);
% Vspread = 1e-7;

Nt_ = [5 9 17 33 65 129];


MC = 1024;

for n = 1:length(Nt_)
    

for k = 1:length(Vspread)

u = normrnd(0, Vspread(k), [1 N]);
v = normrnd(0, 0, [1 N]);


    
for c = 1:MC
    

    
Nt = Nt_(n);


% clear x; clear y; 
% x(1, :) = x0;
% y(1, :) = y0;
% clear z;   



for i = 2:Nt
    x(i, :) = x(i - 1, :) + u .* dt;
    y(i, :) = y(i - 1, :) + v .* dt;
    
%     figure(101); hold on; plot(x(i, :), y(i, :), '+');
    
    D(i, :) = sqrt(x(i, :).^2 + y(i, :).^2);
    z(i) = sum(A0 .* exp(1j .* 4 * pi / lambda .* D(i, :)));
end

Noise = sum(abs(z).^2)./(Nt .* SNR);
sigma_n = sqrt(Noise);

z_obs = z + sigma_n .* (randn(1, Nt));

ZFFT = 1./sqrt(Nt) .* fftshift(fft(z_obs));

v_amb = lambda/(4 * dt);
% dv = lambda/(4 * dt * (Nt - 1));


if (mod(Nt, 2) == 0)
%     hsnew = hs+1;
    vel_axis = linspace(-Nt/2, Nt/2-1, Nt)./Nt .* 2 .* v_amb;
%     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
else
%     hsnew = hs;
    vel_axis = linspace(-v_amb, v_amb, Nt);
end


dv = vel_axis(2) - vel_axis(1); dv_(n) = dv;

%% Plot Doppler Spectrum
% 
% txt = ['Doppler Spectrum'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
% xl = 'velocity [m/s]';
% 
% f = figure(108); hold on; f.Position = [10 10 1000 1000];
% color = 'g';
% yl =  ['Spectrum [dB]'];
% 
% marker = markers(1);
% plott2(vel_axis, db(abs(ZFFT)), xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);
%%

PT = sum(abs(ZFFT).^2  .* dv);
mu(n, k, c) = 1./PT .* sum(vel_axis .* abs(ZFFT).^2 .* dv);
sigma(n, k, c) = sqrt(sum(1./PT .* (vel_axis - mu(n, k, c)).^2 .* abs(ZFFT).^2 .* dv));

end
end
end

Mu_avg = mean(mu, 3);
Sigma_avg = mean(sigma, 3);

%% Plot 
for n = 1:length(Nt_)
%Mean of Means for MC rotations with the same ground truth
txt = ['Retrieved Mean \mu^{Ret}_{V_r} '];
dtext = [' N_t = ', num2str(Nt_(n))];
xl = '\sigma_{V_r, True} [m.sec^{-1}]';
    
marker = markers(1);
f = figure(104); hold on;   f.Position = [10 10 1000 1000];
color = colors(n).c;
yl =  ['\mu^{Ret}_{V_r} [m.s^{-1}]'];
plott2(Vspread, squeeze(Mu_avg(n, :, 1)), xl, yl, txt, 2, dtext, color, marker)

%Mean of Widths for MC rotations with the same ground truth
txt = ['Retrieved Spectrum Width \sigma^{Ret}_{V_r} '];
dtext = [' N_t = ', num2str(Nt_(n))];
xl = '\sigma_{V_r, True} [m.sec^{-1}]';

f = figure(105); hold on; f.Position = [10 10 1000 1000];
% color = 'k';
yl =  ['\sigma^{Ret}_{V_r} [m.s^{-1}]'];
plott2(Vspread, squeeze(Sigma_avg(n, :, 1)), xl, yl, txt, 2, dtext, color, marker)
end
% std of Mean for MC rotations

% txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
% dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
% xl = '\sigma True [m.sec^{-1}]';
% 
% f = figure(106); hold on; f.Position = [10 10 1000 1000];
% color = 'k';
% yl =  ['Std of \mu retrieved [m.s^{-1}], for ', num2str(MC), ' rotations'];
% plott2(Vspread, squeeze(Mu_re_std(1, 1, 1, m, :)), xl, yl, txt, 2, dtext, color, marker)
% 
% % std of Sigma for MC rotations
% 
% txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
% dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
% xl = '\sigma True [m.sec^{-1}]';
% 
% f = figure(107); hold on; f.Position = [10 10 1000 1000];
% color = 'k';
% yl =  ['Std of \sigma retrieved [m.s^{-1}], for ', num2str(MC), ' rotations'];
% plott2(Vspread, squeeze(Sigma_re_std(1, 1, 1, m, :)), xl, yl, txt, 2, dtext, color, marker)
