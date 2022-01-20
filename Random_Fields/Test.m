clear;
close all;
markers = load('../mono/markers.mat');
markers = markers.markers;


N = 1000;

r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 1.8*pi/180;

r = r0 - dr/2 + dr .* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);

% x0 = 10 + 5 .* rand(1, N);
% % y0 = 10 + 5 .* rand(1, N); 
% y0 = zeros(1, N);

x0 = r .* cos(ph); y0 = r .* sin(ph);
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

x(1, :) = x0;
y(1, :) = y0;

% v = linspace(1, 5, N);
u = normrnd(1, 1, [1 N]);
v = normrnd(0, 0, [1 N]);
Nt = 128;

A0 = ones(1, N);
D0 = 0;

lambda = 0.03; dt = 1e-3;

z(1) = sum(A0 .* exp(1j .* 4 * pi / lambda * x0));

% z(1) = 0;
SNR_db = 30; 
SNR = 10^(SNR_db/10);

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


dv = vel_axis(2) - vel_axis(1);

%% Plot Doppler Spectrum

txt = ['Doppler Spectrum'];
dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
xl = 'velocity [m/s]';

f = figure(108); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['Spectrum [dB]'];

marker = markers(1);
plott2(vel_axis, db(abs(ZFFT)), xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);
%%

PT = sum(abs(ZFFT).^2  .* dv);
mu = 1./PT .* sum(vel_axis .* abs(ZFFT).^2 .* dv);
sigma = sqrt(sum(1./PT .* (vel_axis - mu).^2 .* abs(ZFFT).^2 .* dv));
