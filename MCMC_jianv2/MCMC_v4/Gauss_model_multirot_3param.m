clear;
close all;

lambda = 3e-2;
dt = 1e-3;

L = 4 * pi/lambda * dt;

v_amb = lambda/(2 * dt);

% Nu = 20000000;
% 
% U = linspace(0, v_amb, Nu);

Mu = 7;
Sigma = 1;

Gamma = L .* Sigma;
epsilon = 1e-4;
N_min  = 1./Gamma .* ( sqrt(-log(epsilon) - log(1 - exp(-Gamma.^2))) - 1);

Nt = 128;

%K = 0:1:Nt-1;

Nscan = 10000;
M = 1;

K = 100:1:Nt-1;

Mu_width = Mu/20;
Sigma_width = Sigma./20;

Muvec = Mu - Mu_width/2 + Mu_width .* rand(1, M);
Sigmavec = Sigma - Sigma_width/2 + Sigma_width .* rand(1, M);

Ntot = length(K) .* M;


%% Original model

x = [];

%% Generate signal with IFFT

Nu = 200000;

for i = 1:M
    
    U = normrnd(Muvec(i), Sigmavec(i), [1 Nu]);

    for m = 1:length(K)

        y(m) = [sum(exp(1j .* L .* K(m) .* U))];

    end

    x = [x y];

end

%%
% x = [];
% 
% %% Generate signal with IFFT
% 
% U = linspace(0, v_amb, Nt);
% 
% for i = 1:M
% 
% X = 1/sqrt(2 * pi * Sigmavec(i).^2 .* L.^2) .* exp(-(U - Muvec(i)).^2./(2 .* Sigmavec(i).^2));
% 
% x = [x (ifft((X)))];
% 
% end
%%

% x = [];
% for i = 1:M
% 
% % K = [K (i).*Nscan+(0:1:Nt-1)]; % 0+Nscan+1:1:Nscan+Nt];
% 
% x = [x exp(-1/2 * Sigmavec(i).^2 .* L.^2 .* K.^2) .* exp(1j .* L .* K .* Muvec(i))]; %/sqrt(2*pi)];
% 
% end

SNR_db = 30;

SNR = 10^(SNR_db/10);

Noise = sum(abs(x).^2)./(Ntot .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise);

z = x + sigma_n/sqrt(2) .* (randn(1, Ntot) + 1j .* randn(1, Ntot)); % Adding complex noise 


figure; plot(real(z)); hold on; plot(imag(z));

Nfft = 1024;

FFT = fft(z, Nfft);

U = linspace(0, v_amb, Nfft);

figure; plot(U, db(abs(FFT)));

[~, M1, S1] = Mom(FFT, U);

%% MCMC

MU_0 = 1;
SIGMA_0 = 3;
NU_0 = 3000;

Iter = 100000;

% x_mc = [];

MU(1) = MU_0; SIGMA(1) = SIGMA_0; NU(1) = NU_0;
 
for l = 2:Iter

% K = [K (i).*Nscan+(0:1:Nt-1)]; % 0+Nscan+1:1:Nscan+Nt];
Muold = MU(l - 1);
% Munew = normrnd(MU(l - 1), MU(l - 1)/20, 1);
Munew = v_amb .* rand();

Sigmaold = SIGMA(l - 1);
% Sigmanew = normrnd(SIGMA(l - 1), SIGMA(l - 1)./20, 1);
Sigmanew = v_amb/2 .* rand();

Nuold = NU(l - 1);
% Nunew = poissrnd(NU(l - 1));
Nunew = 10000 * rand();

[LLold] = LLGaussv2([Muold Sigmaold Nuold], K, L, M, sigma_n, z);
[LLnew] = LLGaussv2([Munew Sigmanew Nunew], K, L, M, sigma_n, z);

if acceptance(LLold, LLnew)
    MU(l) = Munew;
    SIGMA(l) = Sigmanew;
    NU(l) = Nunew;
else
    MU(l) = Muold;
    SIGMA(l) = Sigmaold;
    NU(l) = Nuold;
end

MUALL(l - 1) = Munew;
SIGMAALL(l - 1) = Sigmanew;
NUALL(l - 1) = Nunew;

end


figure; 
plot(MUALL, 'y', 'LineWidth', 1.2); 
hold on;
plot(MU, 'b', 'LineWidth', 2); 
hold on; 
plot(ones(1, Iter).*Mu, 'r-.', 'LineWidth', 2);
legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;
title({'\mu'}); 

figure; 
plot(SIGMAALL, 'y', 'LineWidth', 1.2); 
hold on;
plot(SIGMA, 'b', 'LineWidth', 2); 
hold on; 
plot(ones(1, Iter).*Sigma, 'r-.', 'LineWidth', 2);
legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;
title({'\sigma'}); 

figure; 
plot(NUALL, 'y', 'LineWidth', 1.2); 
hold on;
plot(NU, 'b', 'LineWidth', 2); 
hold on; 
plot(ones(1, Iter).*Nu, 'r-.', 'LineWidth', 2);
legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;
title({'Number of scatterers N_u'}); 


%% Histogram

Burnin = 100;

figure; histogram(MU(Burnin+1:end)); title({'\mu'}); 

figure; histogram(SIGMA(Burnin+1:end)); title({'\sigma'});

figure; histogram(NU(Burnin+1:end)); title({'Number of scatterers N_u'}); 
