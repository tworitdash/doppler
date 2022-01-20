function [s, v_amb, dv, SNR, hs, N, t, mu, phi_axis, Nphi, sec] = TDG(r, Omega_rpm, n_rot, lambda, BW_deg, PRT, Vxm, Vxs, Vym, Vys, phi_0_deg, theta, n_sig)

Nr = length(r);
BW = BW_deg * pi/180;
Omega = Omega_rpm * 2*pi/60; % rotation speed in rpm


sec = ceil(2*pi/BW);
phi_0 = phi_0_deg * pi/180;


ph = phi_0:BW:(sec)*BW+phi_0;

T = (ph(end)-phi_0)/Omega;


Nf = floor(T/PRT);

hs = floor(Nf/sec);

N = hs*sec;

t = linspace(0, (Nf-1)*PRT, N);

% beta_wind = beta_wind_deg .* pi/180;

v_amb = lambda/(4 * PRT);
dv = lambda/(4*hs*PRT);

s = zeros(n_rot, Nr, N);

phi_axis = mean([ph(1:end-1); ph(2:end)]);
Nphi = length(phi_axis);

for n = 1:n_rot
    for i = 1:Nr
            Vx = normrnd(Vxm(i), Vxs(i), [1 N]);
            Vy = normrnd(Vym(i), Vys(i), [1 N]);
    %         mu = V .* cos((beta_wind(i) - Omega .* t - phi_0)) .* cos(theta);
            mu(n, i, :) = Vx .* cos(Omega.* t - phi_0) .* cos(theta) + Vy .* sin(Omega.*t - phi_0) .* cos(theta);

            [s(n, i, :), SNR(n, i)] = TD_generatorV2(mu(n, i, :), lambda, t, n_sig);
    end
end
