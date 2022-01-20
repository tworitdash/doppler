function [data_f, data, X, Theta, vel_axis] = DS(SNR, m0, mu, sigma, n, v_amb)

%% N is the number of Doppler bins permitted for a specific rotation speed and radar beam width 
if mod(n, 2) == 0
    
    axis = linspace(-n/2, n/2-1, n)/(n);
    vel_axis = 2 * v_amb * axis;
    
else
    
    vel_axis = linspace(-v_amb, v_amb, n);
    
end

X = rand(1, n);
Theta = 2 .* pi * rand(1, n);

if sigma < 0.02
    [~, idx1] = (min(abs(vel_axis - mu)));
    S_ = dirac(vel_axis - vel_axis(idx1)); 
    idx = S_ == Inf;
    S_(idx) = 1;
else
S_ = m0/sqrt(2*pi*sigma^2) * exp(-(vel_axis - mu).^2/(2*sigma^2));
end

Noise_full = sum(S_) ./ (n .* SNR); % Noise Power

P_full = -(S_ + Noise_full) .* log(X); % Power spectrum 

data_f = sqrt(P_full) .* exp(1j .* Theta);
data = ifft(fftshift(sqrt(n) .* data_f));

end