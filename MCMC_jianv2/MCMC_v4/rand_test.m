clear;
mu = 2;
sig = 1;

% f = @(x) 1/sqrt(2*pi*sig^2) .* exp(-(x - mu).^2./(2*sig.^2));
f = @(x) erfinv(2*x-1)*sig*sqrt(2) +mu;

MIN = mu - 3*sig;
MAX = mu + 3*sig;

N = 100000;

u = MIN + (MAX - MIN) .* rand(1, N);
figure; histogram(u);

g = f(u);

figure; histogram(g);