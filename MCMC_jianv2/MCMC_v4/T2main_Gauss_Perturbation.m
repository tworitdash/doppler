clear;

Nu = 10000;
Mu = pi; 
% Sigma = 0.5 * 4 * pi / 0.03 * 1e-3;
Sigma = pi/10;
U_ = normrnd(Mu, Sigma, [1 Nu]);

Psi_ = -pi + 2 * pi * rand([1 Nu]);
% Psi_ = 0 .* rand([1 Nu]);

[Psim, Psin] = meshgrid(Psi_, Psi_);

[U1, U2] = meshgrid(U_, U_);

Np = 16;
Ng = 0;
Ns = 20;
K_ = [];

for i = 1:Ns
    K_ = [K_ (i-1)*(Np + Ng):(i-1)*(Np+Ng)+Np-1];
end

Nt = length(K_);

for ki = 1:Nt
    k = K_(ki);
    x(ki) = sum( exp(1j .* U_ .* k + 1j .* Psi_) );
end

SNR_db = 60;

SNR = 10^(SNR_db/10);

Noise = sum(abs(x).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise);

z = x + sigma_n/sqrt(2) .* (randn(1, Nt) + 1j .* randn(1, Nt)); % Adding complex noise 


u = linspace(0, 2 * pi, Np);

for i = 1:Ns
    Z(i, :) = fft(z((i-1)*(Ng+Np)+1:(i-1)*(Ng+Np)+Np));
    figure(1); hold on; plot(u, (abs((Z(i, :)))));
end


for ui = 1:Np
%     figure(1000+ui); histogram(real(Z(:, ui))); title([num2str(ui)]);
    mu(ui) = mean(real(Z(:, ui)));
    sig(ui) = std(real(Z(:, ui)));
end

figure(2); hold on; plot(u, mu); 

figure(3); hold on; plot(u, sig);
