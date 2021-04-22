clear; close all;

f = 10e9;
n = 1024;

t = linspace(eps, 1/f, n);

x = exp(1j .* 2 .* pi .* f .* t);

x_f = fft(fftshift(x))/sqrt(n);

x_t = ifft(fftshift(x_f .* sqrt(n)));

x_t_modified = abs(x_t) .* exp(1j .* angle(x_t) .* (-1));
x_f_modified = fft(fftshift(x_t_modified))/sqrt(n);

f_ = -10e9;
n_ = n;

t_ = linspace(eps, 1/f_, n_);

x_ = exp(1j .* 2 .* pi .* f_ .* t_);

x_f_ = fft(fftshift(x_))/sqrt(n);



figure; plot(db(abs(x_f_modified))); hold on; plot(db(abs(x_f))); hold on; plot(db(abs(x_f_)));