clear; 
close all;

lambda = 3e-2;
dt = 1e-3;
L = 2/lambda;

m0 = 1;

mu = 7.5*L;

f_amb = 1/dt;

N = 2048;

sigma = 1*L;

SNR_db = 30;
SNR = 10^(SNR_db/20);


[data, vel, data_nf, Noise, ~] = DS(SNR, m0, mu, sigma, N, f_amb);

FFT = abs(fftshift(fft(data)));

figure; plot(vel*lambda/2, db(FFT));

[PT, M, S] = Mom(FFT, vel*lambda/2);


figure; plot(real(data)); hold on; plot(imag(data));
figure; plot(real(data_nf)); hold on; plot(imag(data_nf));

%% Known quantities

sigma_n = sqrt(Noise);

%% MCMC

MU_0 = 7.5*L;
SIGMA_0 = 1*L;

Iter = 1000;

% x_mc = [];

MU(1) = MU_0; SIGMA(1) = SIGMA_0;
 
for l = 2:Iter

% K = [K (i).*Nscan+(0:1:Nt-1)]; % 0+Nscan+1:1:Nscan+Nt];
Muold = MU(l - 1);
% Munew = normrnd(MU(l - 1), MU(l - 1)/20, 1);
Munew = f_amb .* rand();

Sigmaold = SIGMA(l - 1);
% Sigmanew = normrnd(SIGMA(l - 1), SIGMA(l - 1)./20, 1);
Sigmanew = f_amb/2 .* rand();

[LLold] = LLGaussifft(m0, Muold, Sigmaold, sigma_n, data, N, f_amb, SNR);
[LLnew] = LLGaussifft(m0, Munew, Sigmanew, sigma_n, data, N, f_amb, SNR);

if acceptance(LLold, LLnew)
    MU(l) = Munew;
    SIGMA(l) = Sigmanew;
    
    
else
    MU(l) = Muold;
    SIGMA(l) = Sigmaold;
end

MUALL(l - 1) = Munew;
SIGMAALL(l - 1) = Sigmanew;


end


figure; 
plot(MUALL, 'y', 'LineWidth', 1.2); 
hold on;
plot(MU, 'b', 'LineWidth', 2); 
hold on; 
plot(ones(1, Iter).*mu, 'r-.', 'LineWidth', 2);
legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;

figure; 
plot(SIGMAALL, 'y', 'LineWidth', 1.2); 
hold on;
plot(SIGMA, 'b', 'LineWidth', 2); 
hold on; 
plot(ones(1, Iter).*sigma, 'r-.', 'LineWidth', 2);
legend({'MCMC All samples', 'MCMC', 'Ground Truth'});  grid on;

%% Histogram

Burnin = 100;

figure; histogram(MU(Burnin+1:end)); 

figure; histogram(SIGMA(Burnin+1:end));

