clear;
% close all;

Sigma = linspace(0.1, 4, 1000); 

Ntvec = 1:1:1024;

for s = 1:length(Sigma)

sigma = Sigma(s);

lambda = 3e-2; 
dt = 1e-3;

L = 4 * pi * dt/lambda;

G(s) = sqrt(L.^2 .* sigma.^2);

for i = 1:length(Ntvec)

Nt = Ntvec(i);


k = 0:1:Nt-1;

S(i) = sum(exp(-k.^2.*L.^2.*sigma.^2));

end

for m = 2:length(Ntvec)
    dSdN(m-1) = (S(m) - S(m-1));
end

H = (dSdN < 1e-5);

K_conv = k(H == 1);

Nt_conv(s) = k(K_conv(1));

end

figure; plot(Sigma, Nt_conv); grid on; hold on; plot(Sigma, pi./G)

% figure(1); hold on; plot(Ntvec, S, 'DisplayName', ['L^2\sigma^2 = ', num2str(G^2)]); legend;

