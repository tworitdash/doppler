function [ZFFT, PT, mu, sigma, vel_axis, dv, v_amb] = Spec2(z_model, Nt, dt, lambda, SNR_db, isplot, mi, ci)
markers = load('markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;


ZFFT = 1./sqrt(Nt) .* fftshift(fft(z_model));
v_amb = lambda/(4 * dt);

% v_amb = 1/(2 * dt);
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

PT = sum(abs(ZFFT).^2  .* dv);
mu = 1./PT .* sum(vel_axis .* abs(ZFFT).^2 .* dv);
sigma = sqrt(sum(1./PT .* (vel_axis - mu).^2 .* abs(ZFFT).^2 .* dv));

%% Plot Doppler Spectrum
if isplot
txt = ['Doppler Spectrum'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]'];
xl = 'velocity [m/s]';

f = figure(108); hold on; f.Position = [10 10 1000 1000];
color = colors(ci).c;
yl =  ['Spectrum [dB]'];

marker = markers(mi);
plott2(vel_axis, db(abs(ZFFT)), xl, yl, txt, 2, dtext, color, marker); % hold on; histogram(v);
end
%%



end