clear;
close all;

L = pi/4;
L1 = pi/4;

N = 1;
% 
% F = -L + 2 * L * rand([1 N]);
% phi = -pi/4 + pi/2 * rand([1 N]);

% figure; histogram(F+phi, 100);

for k = 1:2^20
    F = -L + 2 * L * rand([1 N]);
    phi = -L1 + 2* L1 * rand([1 N]);

%     figure(100); hold on; histogram(F+phi, 100);
    Z(k) = sum(exp(1j .* (F + phi)));
end

figure; histogram(real(Z)); hold on; histogram(imag(Z));

mr = mean(real(Z)); sr = std(real(Z));
mi = mean(imag(Z)); si = std(imag(Z));

Erz = N./(2*L) .* ( (2 * (4 - pi) * sin(L) + pi * cos(L))./(pi * sqrt(2)) + 2 * sin(L - L1) ); 