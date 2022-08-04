clear;


lambda = 3e-2;
dt = 1e-3;

L = 4 * pi/lambda * dt;

v_amb = lambda/(2 * dt);

Nu = 20000000;

U = linspace(0, v_amb, Nu);

Mu = 7.5;
Sigma = 1.2;


Spectrum = 1./sqrt(2 * pi * Sigma.^2 .* L.^2) .* exp(-(U - Mu).^2./(2 .* Sigma.^2));

% figure; plot(U, db(Spectrum));


Nt = 128;

u_ = linspace(0, v_amb, Nt);

for i = 1:length(u_)
    [u(i), indx(i)] = q1(U, u_(i), U);
end

Spectrum_Nt = Spectrum(indx);

V = normrnd(Mu, Sigma, [1 Nu]);

% figure; histogram(V);

for k = 0:Nt-1

 % sum( Spectrum_Nt .* exp(-1j .* L .* u .* k));

IFFT_(k+1) = sum( exp(1j .* L .* V .* k) ) ./ sum(Spectrum);

end

% IFFT = ifft(fftshift(Spectrum_Nt));

figure; plot(real(IFFT_)./max(real(IFFT_))); hold on; plot(imag(IFFT_)./max(imag(IFFT_)));


IFFT_a = exp(-1./2 .* Sigma.^2 .* L.^2 .* (0:1:Nt-1).^2) .* exp(  1j .* Mu .* L .* (0:1:Nt-1))./sqrt(2 * pi); 
% figure; plot(u, abs((fft(IFFT_a)))); 


hold on; plot(real(IFFT_a)./max(real(IFFT_a))); hold on; plot(imag(IFFT_a)./max(imag(IFFT_a)));

Nfft = 128;
Ufft = linspace(0, v_amb, Nfft);

FFT_a = abs((fft(IFFT_a, Nfft)));
FFT_ = abs((fft(IFFT_, Nfft)));

figure; plot(Ufft, db((FFT_ ./ max(FFT_)))); 

hold on; plot(Ufft, db((FFT_a ./ max(FFT_a)))); 

[~, M1, S1] = Mom(FFT_, Ufft);
[~, M2, S2] = Mom(FFT_a, Ufft);

