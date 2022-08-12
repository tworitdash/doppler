clear;
% close all;

lambda = 3e-2;
dt = 1e-3;

L = 4 * pi/lambda * dt;

v_amb = lambda/(2 * dt);

% Nu = 20000000;
% 
% U = linspace(0, v_amb, Nu);

Mu = 7.5;
Sigma = 1;

Gamma = L .* Sigma;
epsilon = 1e-4;
N_min  = 1./Gamma .* ( sqrt(-log(epsilon) - log(1 - exp(-Gamma.^2))) - 1);

Nt = 128;

%K = 0:1:Nt-1;

Nscan = 100;
M = 1;
K = [];
for n = 1:M
    K = [K Nscan*(n-1):1:Nscan*(n-1)+Nt-1];
end

Mu_width = Mu/20;
Sigma_width = Sigma./20;

Muvec = Mu; % - Mu_width/2 + Mu_width .* rand(1, M);
Sigmavec = Sigma; % - Sigma_width/2 + Sigma_width .* rand(1, M);

Ntot = Nt .* M;


%% Original model

x = [];

%% Generate signal with IFFT

Nu = 100000;
Phi = 2 * pi .* rand([1 Nu]);
% Phi = 0;

for i = 1:M
    
    U = normrnd(Muvec(i), Sigmavec(i), [1 Nu]);
%     U = Muvec(i);

    for m = 1:Nt

        x((i-1).*Nt+m) = [sum(exp(1j .* L .* K((i-1).*Nt+m) .* U + 1j .* Phi) )];
%         y(m) = [sum( exp(1j .* Phi) )];

    end

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
% x = [x exp(-1/2 * Sigmavec(i).^2 .* L.^2 .* K.^2) .* exp(1j .* L .* K .* Muvec(i))]; %/sqrt(2*pi)];
% 
% end

SNR_db = 30;

SNR = 10^(SNR_db/10);

Noise = sum(abs(x).^2)./(Ntot .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise);

z = x + sigma_n/sqrt(2) .* (randn(1, Ntot) + 1j .* randn(1, Ntot)); % Adding complex noise 


figure; plot(real(x)); hold on; plot(imag(x)); grid on;

Nfft = 128;

FFT = fft(z, Nfft);

U = linspace(0, v_amb, Nfft);

figure; plot(U, db(abs(FFT)));

[~, M1, S1] = Mom(FFT, U);


%% LL Test

MU_0 = Mu;

Iter = 400;
MC = 2048;

SIGMA_ = [linspace(0.8, 1, 100) linspace(1, 4, 10000)] ;

% x_mc = [];

MU(1) = MU_0; 

parfor i = 1:length(SIGMA_)
disp(i);
SIGMA = SIGMA_(i);
 
[LL(i)] = LLGaussrandom(MU_0, SIGMA, K, L, M, sigma_n, z, Nu, Nt, MC);

end

figure; plot(SIGMA_, LL);