f = 1e3/2;

T = 2/f;


N = 1025;


t = linspace(eps, T, N);

x = sin(2*pi*f*t);

% plot(t, x);

y = fft(x, N);

x1 = ifft(y, N);

y1 = fft(x1, N);

figure; plot(abs(y).^2); hold on; plot(abs(y1).^2);