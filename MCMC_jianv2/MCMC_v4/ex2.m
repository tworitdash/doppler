clear;
% close all;

sigma = 0.2; 

Ntvec = 0:1:512;

lambda = 3e-2; 
dt = 1e-4;

L = 4 * pi * dt/lambda;

G = sqrt(L.^2 .* sigma.^2);

for i = 1:length(Ntvec)

Nt = Ntvec(i);


k = 0:1:Nt-1;

S(i) = sum(exp(-k.^2.*L.^2.*sigma.^2));

end


figure(1); hold on; plot(Ntvec, S, 'DisplayName', ['L^2\sigma^2 = ', num2str(G^2)]); legend;
