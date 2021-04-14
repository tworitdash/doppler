function [data, data_f] = DS_simulator(SNR, m0, mu, sigma, n, v_amb)


vel_axis = linspace(-v_amb, v_amb, n);

X = rand(1, n);
Theta = 2 .* pi * rand(1, n);

if sigma < 0.02
    [~, idx1] = min(abs(vel_axis - mu));
    S = dirac(vel_axis - vel_axis(idx1)); 
    idx = S == Inf;
    S(idx) = 1;
else
    S = m0/sqrt(2*pi*sigma^2) * exp(-(vel_axis - mu).^2/(2*sigma^2));
end

% n = length(S); % Length of number of frequencies

N = sum(S) ./ (n .* SNR); % Noise Power

P = -(S + N) .* log(X); % Power spectrum 
data_f = sqrt(P);

% figure(1); hold on; plot(vel_axis, sqrt(abs(P))); grid on;

data = ifft(fftshift(sqrt(n) .* sqrt(P) .* exp(1j .* Theta)));% complex time domain signal 
% data = ifft(fftshift(sqrt(P) .* exp(1j .* Theta)));% complex time domain signal 

% figure; plot(abs(data)); 

end