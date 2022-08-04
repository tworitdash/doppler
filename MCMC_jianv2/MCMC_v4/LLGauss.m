
function LL = LLGauss(MU, SIGMA, K, L, M, sigma_n, z, Nu)
x_mc = [];

Phi = -pi + 2 .* pi .* rand([1 Nu]);

for i = 1:M
    
x_mc = [x_mc Nu .* MU .* exp(-1/2 * SIGMA.^2 .* L.^2 .* K.^2) .* exp(1j .* L .* K .* MU)./(2 * pi)];

end

% figure; plot(real(x_mc)); hold on; plot(imag(x_mc));

X = [real(x_mc) imag(x_mc)];
Z = [real(z) imag(z)];

LL = -[Z - X] * [Z - X].' / (2 * sigma_n.^2);

end