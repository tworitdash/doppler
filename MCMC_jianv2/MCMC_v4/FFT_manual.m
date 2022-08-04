clear;
markers = load('markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

dt = 1e-3;
lambda = 3e-2;

L = 4 * pi / lambda * dt;

Nt = 5;

Mu = 7.5;
k = 0:1:Nt-1;

u = linspace(0, lambda/(2 * dt), 128);
Sigmavec = linspace(eps, 4, 100);

for s = 1:length(Sigmavec)

Sigma = Sigmavec(s);


for i = 1:length(u)

FFT(i) = sqrt(25600) .* sum(exp(-Sigma.^2 .* L.^2 .* k.^2) .* exp(1j .* k .* (L * Mu - L .* u(i))));

end

% figure(108); plot(u, db((abs(FFT)).^2)-3);

dv = u(2) - u(1);

dvres = lambda/(2 * dt)/(Nt - 1);

PT(s) = sum(abs(FFT).^2  .* dv);
mu(s) = 1./PT(s) .* sum(u .* abs(FFT).^2 .* dv);
sigma(s) = sqrt(sum(1./PT(s) .* (u - mu(s)).^2 .* abs(FFT).^2 .* dv));


% 
% PT = sum(abs(ZFFT).^2  .* dv);
% mu = 1./PT .* sum(vel_axis .* abs(ZFFT).^2 .* dv);
% sigma = sqrt(sum(1./PT .* (vel_axis - mu).^2 .* abs(ZFFT).^2 .* dv));
end

txt = ['\sigma_{True} vs \sigma_{re}'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]'];

dtext = ['Nt = ', num2str(Nt)];
% ', dv = ', num2str(dv_obs), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]', ...
%    ', \mu_{True} = ', num2str(mu_true), ' [m/s]', ', \sigma_{True} = ', num2str(sigma_true), ' [m/s]'];

xl = '\sigma_{True} [m/s]';

f = figure(1); hold on; f.Position = [10 10 1000 1000];
color = colors(4).c;
yl =  ['\sigma_{re}'];

marker = markers(1);
plott2(Sigmavec, sigma, xl, yl, txt, 2, dtext, color, marker);

marker = markers(2);

plott2(pi./(L .* Nt) .* ones(1, length(Sigmavec)), ...
    sigma, xl, yl, txt, 2, dtext, color, marker);