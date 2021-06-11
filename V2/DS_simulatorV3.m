function [data, data_f, X, Theta] = DS_simulatorV3(SNR, m0, mu, sigma, N, v_amb)

%% N is the number of Doppler bins permitted for a specific rotation speed and radar beam width 
if mod(N, 2) ~= 0
    axis = linspace(-(N - 1)/2, (N - 1)/2, N)/(N - 1);
    vel_axis = 2 * v_amb * axis;
else
    axis = linspace(-N/2, N/2 - 1, N)/N;
    vel_axis = 2 * v_amb * axis;
end
% vel_axis = linspace(-v_amb, v_amb, n);
dv = vel_axis(2) - vel_axis(1);

X = rand(1, N);
Theta = 2 .* pi * rand(1, N);


if sigma < 0.02
    [~, idx1] = (min(abs(vel_axis - mu)));
    S_ = dirac(vel_axis - vel_axis(idx1)); 
    idx = S_ == Inf;
    S_(idx) = 1;
else
    S_ = m0/sqrt(2*pi*sigma^2) * exp(-(vel_axis - mu).^2/(2*sigma^2));
end

Noise = sum(S_) ./ (N .* SNR); % Noise Power


P = -(S_ + Noise) .* log(X); % Power spectrum 
data_f = sqrt(P) .* exp(1j .* Theta);
data = ifft(fftshift(sqrt(N) .* data_f));

end