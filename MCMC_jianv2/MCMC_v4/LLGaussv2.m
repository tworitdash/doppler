function LL = LLGaussv2(x, K, L, M, sigma_n, z)
x_mc = [];

MU = x(1);
SIGMA = x(2);
NU = x(3);
ep = 1e-3;

for i = 1:M
    
x_mc = [x_mc NU.^2 .* MU .* exp(-1/2 * SIGMA.^2 .* L.^2 .* K.^2) .* exp(1j .* L .* K .* MU)./(2 * pi) .* ep];

end



X = [real(x_mc) imag(x_mc)];
Z = [real(z) imag(z)];

LL = -[Z - X] * [Z - X].' / (2 * sigma_n.^2);

end