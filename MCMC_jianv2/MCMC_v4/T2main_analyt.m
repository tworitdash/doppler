clear;
% close all;
Nu = 10000;

Mu = pi; 
Sigma = 1 * 4 * pi / 0.03 * 1e-3;
U_ = normrnd(Mu, Sigma, [1 Nu]);

Np = 5;
Ng = 95;
Ns = 2;
K_ = [];
for i = 1:Ns
    K_ = [K_ (i-1)*(Np + Ng):(i-1)*(Np+Ng)+Np-1];
end


[U, K] = meshgrid(U_, K_);

Psi = -pi + 2 * pi * rand([1, Nu]);
N = 10000;
u = linspace(0, 2*pi, N);

for ui = 1:N
    LL(ui) = - sum( (sum(cos(U_.*K + Psi), 2).' - cos(u(ui) .* K_)).^2 + (sum(sin(U_.*K + Psi), 2).' - sin(u(ui) .* K_)).^2 );
end


figure(1); hold on; plot(u, LL - max(LL));