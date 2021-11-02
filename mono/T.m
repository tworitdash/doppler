%% Radar parameters
clear;
lambda = 0.03;

hs = 5;
sec = 200;
N = hs*sec;

v = linspace(-7.5, 7.5, hs);
phi = linspace(0, 360, N) .* pi/180;
phi_axis = linspace(0, 360, sec);

mu = 5 * cos(phi);
sigma = 0.1 * cos(phi);

SNR = 10^(30/10);
m0 = 1;
figure;
for i = 1:sec
    idx = (i - 1)*hs+1:i*hs;
    [S] = Gauss_PDF(m0, mu(idx), sigma(idx), v);
    X = rand(1, hs);
    Theta = 2 * pi * rand(1, hs);
    N = sum(S)./(hs * SNR);
    P = -(S + N).*log(1 - X);
    sf = sqrt(P) .* exp(1j .* Theta);
    hold on; plot(v, abs(sf));
    s(idx) = ifft(ifftshift(sqrt(hs) .* sf));
end 

Sre = spectrogram(s, hs, 0, hs, N/60);

figure; surface(v, phi_axis, db(abs(Sre)).'); shading flat; 