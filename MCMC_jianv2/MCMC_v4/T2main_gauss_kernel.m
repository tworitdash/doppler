clear;

Nu = 10000;
Mu = pi; 
% Sigma = 0.5 * 4 * pi / 0.03 * 1e-3;
Sigma = 0.3;
U_ = normrnd(Mu, Sigma, [1 Nu]);

% Psi_ = -pi + 2 * pi * rand([1 Nu]);
Psi_ = 0 .* rand([1 Nu]);

[Psim, Psin] = meshgrid(Psi_, Psi_);

[U1, U2] = meshgrid(U_, U_);

Np = 128;
Ng = 0;
Ns = 2;
K_ = [];

for i = 1:Ns
    K_ = [K_ (i-1)*(Np + Ng):(i-1)*(Np+Ng)+Np-1];
end


S = linspace(0, 2*pi, 100);
M = linspace(0, 2 * pi, 100);

for s = 1:length(S)
    for m = 1:length(M)

            for ki = 1:length(K_)
                        k = K_(ki);
    
                        Lk(m, s, ki) = sum( exp(1j .* ((U_ - M(m)) .* k + Psi_) - S(s).^2 .* k.^2 / 2) );
            end
        L(s, m) = sum(Lk(m, s, :), 3);
    end
end


figure; surface(M, S, db(abs(L))); shading flat;
figure; plot(S, abs(L(:, 51)))
figure; plot(M, db(abs(L(1, :))));