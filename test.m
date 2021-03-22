f = -10:0.1:10;

omega = 2 * pi * f;
a = 1;

sigma_ep = 1;
sigma_z = 1;

F_omega = (sigma_ep.^2)./(1 - a .* exp(-1j .* omega)).^2;

figure; plot(f, abs(F_omega));