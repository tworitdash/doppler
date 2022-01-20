%% 1D experiment
clear;

v = 5;
dt = 1e-3;
dx = v*dt;
sigma_l = 3;

x = eps:dx:10;

A(1, :) = exp(-(x - x(1)).^2./(2 .* sigma_l.^2));

Nt = 1024;

t = eps:dt:Nt*dt;
% figure(1); plot(x, A(1, :)); 

D(1, :) = 0;

for i = 2:Nt
    A(i, :) = exp(-(x - v.*i.*dt).^2./(2 .* sigma_l.^2));
    D(i, :) = D(i - 1, :) + v .* dt;
end

lambda = 3e-2;

for i = 1:Nt
    z(i) = sum(A(i, :) .* exp(1j .* 4*pi/lambda .* D(i, :)));
end

ZFFT = 1./sqrt(Nt) .* fftshift(fft(z));

v_amb = lambda/(4 * dt);
vel_axis = linspace(-v_amb, v_amb, Nt);

figure; plot(vel_axis, abs(ZFFT).^2);

%% Experiment 2:
clear;

v = 5;
dt = 1e-3;
dx = v*dt;
x0 = eps:dx:10;
lambda = 3e-2;

Nt = 1024;

A = ones(size(x0));

x(1, :) = x0;
for i = 2:Nt
    x(i, :) = x(i - 1, :) + v .* dt;
    z(i) = sum(A .* exp(1j .* 4*pi/lambda .* x(i, :)));
end

ZFFT = 1./sqrt(Nt) .* fftshift(fft(z));

v_amb = lambda/(4 * dt);
vel_axis = linspace(-v_amb, v_amb, Nt);

figure; plot(vel_axis, abs(ZFFT).^2);
