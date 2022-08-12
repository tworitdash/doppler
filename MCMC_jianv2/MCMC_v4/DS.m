function [data, vel_axis, data_f_Sig, Noise, P] = DS(SNR, m0, mu, sigma, N, f_amb)

vel_axis = linspace(0, f_amb, N);

X = rand(1, N);
Theta = 2 .* pi * rand(1, N);

for i = 1:N
    S_(i) = m0./sqrt(2*pi*(sigma).^2) .* exp(-(vel_axis(i) - mu).^2./(2*(sigma).^2));
end

% figure; plot(abs(S_))
% end
%% Doviak

Noise = sum(S_) ./ (N .* SNR); % Noise Power

P_wNoise = -(S_ + Noise) .* log(1 - X); % Power spectrum 

data_f = sqrt(P_wNoise) .* exp(1j .* Theta);

data = ifft(fftshift(sqrt(N) .* data_f));

% figure; plot(db(abs(data_f)));

%% Only signal and no noise

P = -(S_); % Power spectrum 

data_f_Sig = ifft(ifftshift(sqrt(N) .* sqrt(P) .* log(1 - X) .* exp(1j .* Theta) )); % .* exp(1j .* Theta);
 


end