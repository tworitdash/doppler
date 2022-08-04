clear;
close all;

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

Nt = 5;

%K = 0:1:Nt-1;

Nscan = 1000;
M = 1;

K = 0:1:Nt-1;

Mu_width = Mu/20;
Sigma_width = Sigma./20;

Muvec = Mu - Mu_width/2 + Mu_width .* rand(1, M);
Sigmavec = Sigma - Sigma_width/2 + Sigma_width .* rand(1, M);

Ntot = Nt .* M;


%% Original model

x = [];

%% Generate signal with IFFT

Nu = 1000;

for i = 1:M
    
    U = normrnd(Muvec(i), Sigmavec(i), [1 Nu]);

    for m = 1:Nt

        y(m) = [sum(exp(1j .* L .* K(m) .* U))]./(sum(U));

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

Nfft = 128;

FFT = fft(z, Nfft);

U = linspace(0, v_amb, Nfft);

figure; plot(U, db(abs(FFT)));

[~, M1, S1] = Mom(FFT, U);



%% Likelihood test

Mull = linspace(0, v_amb, 100);

Sigmall = linspace(0, 4, 100);

for p = 1:length(Mull)
    for q = 1:length(Sigmall)
        
    [LL(p, q)] = LLGauss(Mull(p), Sigmall(q), K, L, M, sigma_n, z);
    
    end
end
x = Mull;
y = Sigmall;
z = LL;

xl = [' Mean Velocity \mu [m/s] ' ];
yl = [' Velocity Spectrum Width \sigma [m/s] ' ];
zl = ['Likelihood [dB]' ];
tl = ['Likelihood vs \mu and \sigma' ];
surplot_pcolor(x, y, z.', xl, yl, zl, tl)

[Mumax, Muidx] = max(LL);

[Sigmamax, Sigmaidx] = max(Mumax);

Mumax_ = Mull(Muidx(1)); Sigmamax_ = Sigmall(Sigmaidx); 


% x = input_info.N_pulsevec;
% y = input_info.velocity.u.sigmavec;
% z = Sigma_error;
% 
% xl = [' Number of Pulses ' ];
% yl = [' True velocity width u_{sigma}' ];
% zl = ['Error in velocity width estimation u_{\sigma_{re}}' ];
% tl = ['Error in velocity width estimation u_{\sigma_{re}}' ];
% surplot_pcolor(x, y, z.', xl, yl, zl, tl)

%% MCMC
% 
% MU_0 = 1;
% SIGMA_0 = 3;
% 
% Iter = 10000;
% 
% % x_mc = [];
% 
% MU(1) = MU_0; SIGMA(1) = SIGMA_0;
%  
% for l = 2:Iter
% 
% % K = [K (i).*Nscan+(0:1:Nt-1)]; % 0+Nscan+1:1:Nscan+Nt];
% Muold = MU(l - 1);
% Munew = normrnd(MU(l - 1), MU(l - 1)/20, 1);
% 
% Sigmaold = SIGMA(l - 1);
% Sigmanew = normrnd(SIGMA(l - 1), SIGMA(l - 1)./20, 1);
% 
% [LLold] = LLGauss(Muold, Sigmaold, K, L, M, sigma_n, z);
% [LLnew] = LLGauss(Munew, Sigmanew, K, L, M, sigma_n, z);
% 
% if acceptance(LLold, LLnew)
%     MU(l) = Munew;
%     SIGMA(l) = Sigmanew;
%     
%     
% else
%     MU(l) = Muold;
%     SIGMA(l) = Sigmaold;
% end
% 
% MUALL(l - 1) = Munew;
% SIGMAALL(l - 1) = Sigmanew;
% 
% 
% end
% 
% 
% figure; 
% plot(MUALL, 'y', 'LineWidth', 1.2); 
% hold on;
% plot(MU, 'b', 'LineWidth', 2); 
% hold on; 
% plot(ones(1, Iter).*Mu, 'r-.', 'LineWidth', 2);
% legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;
% 
% figure; 
% plot(SIGMAALL, 'y', 'LineWidth', 1.2); 
% hold on;
% plot(SIGMA, 'b', 'LineWidth', 2); 
% hold on; 
% plot(ones(1, Iter).*Sigma, 'r-.', 'LineWidth', 2);
% legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;
% 
% %% Histogram
% 
% Burnin = 100;
% 
% figure; histogram(MU(Burnin+1:end)); 
% 
% figure; histogram(SIGMA(Burnin+1:end));