function LL = LLGaussrandom(MU, SIGMA, K, L, M, sigma_n, z, Nu, Nt, MC)
x_mc = [];
LL_ = zeros(1, MC);

for mc = 1:MC

Phi = 2 .* pi .* rand([1 Nu]);
% Phi = 0;
U = normrnd(MU, SIGMA, [1 Nu]);
Ntot = Nt * M;

parfor i = 1:Ntot
    
% x_mc = [x_mc Nu .* MU .* exp(-1/2 * SIGMA.^2 .* L.^2 .* K.^2) .* exp(1j .* L .* K .* MU)./(2 * pi)];

%     for m = 1:Nt

        x_mc(i) = [sum(exp(1j .* L .* K(i) .* U + 1j .* Phi) )];
%         y(m) = [sum( exp(1j .* Phi) )];

%     end

end



% figure; plot(real(x_mc)); hold on; plot(imag(x_mc));

X = [real(x_mc) imag(x_mc)];
Z = [real(z) imag(z)];

LL_(mc) = -[Z - X] * [Z - X].' / (2 * sigma_n.^2);

end
% vr = linspace(0, 15, Ntot);
% 
% figure(2); hold on; plot(vr, db(abs(fft(x_mc))));

LL = mean(LL_);
end