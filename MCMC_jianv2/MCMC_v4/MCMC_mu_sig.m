clear;


lambda = 3e-2;
dt = 1e-3;

L = 4 * pi/lambda * dt;

v_amb = lambda/(2 * dt);

Nu = 20000000;

U = linspace(0, v_amb, Nu);

Mu = 7.5;
Sigma = 2;

Nt = 128;

%K = 0:1:Nt-1;

Nscan = 100;
M = 2;

K = 0:1:Nt-1;

for i = 1:M-1

K = [K (i).*Nscan+(0:1:Nt-1)]; % 0+Nscan+1:1:Nscan+Nt];

end

Ntnew = length(K);

x = exp(-1/2 .* Sigma.^2 .* L.^2 .* K.^2) .* exp(1j .* Mu .* L .* K)./(sqrt(2*pi));

SNR_db = 30000;

SNR = 10^(SNR_db/10);

Noise = sum(abs(x).^2)./(Ntnew .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise);

z = x + sigma_n/sqrt(2) .* (randn(1, Ntnew) + 1j .* randn(1, Ntnew)); % Adding complex noise 

Nfft = 128;

FFT = fft(z, Nfft);

u = linspace(0, v_amb, Nfft);

figure(6); hold on; plot(u, db(abs(FFT)));


[~, M1, S1] = Mom(FFT, u);

Gamma = L .* Sigma;
epsilon = 1e-4;
N_min  = 1./Gamma .* ( sqrt(-log(epsilon) - log(1 - exp(-Gamma.^2))) - 1);

