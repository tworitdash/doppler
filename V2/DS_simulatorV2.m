function [data, data_f, data_f_full, data_full, X_full, Theta_full] = DS_simulatorV2(SNR, m0, mu, sigma, n, v_amb, N)

%% N is the number of Doppler bins permitted for a specific rotation speed and radar beam width 


vel_axis = linspace(-v_amb, v_amb, n);
dv = vel_axis(2) - vel_axis(1);

X_full = rand(1, n);
Theta_full = 2 .* pi * rand(1, n);

% if sigma < 0.02
%     [~, idx1] = (min(abs(vel_axis - mu)));
%     S_ = dirac(vel_axis - vel_axis(idx1)); 
%     idx = S_ == Inf;
%     S_(idx) = 1;
% else
S_ = m0/sqrt(2*pi*sigma^2) * exp(-(vel_axis - mu).^2/(2*sigma^2));
% end

Noise_full = sum(S_) ./ (n .* SNR); % Noise Power

P_full = -(S_ + Noise_full) .* log(X_full); % Power spectrum 

data_f_full = sqrt(P_full) .* exp(1j .* Theta_full);
data_full = ifft(fftshift(sqrt(n) .* data_f_full));

vel_axis_permitted = linspace(-v_amb, v_amb, N);

for i = 1:N
    [~, indices(i)] = min(abs(vel_axis - vel_axis_permitted(i)));
end

idx_for_integral = round(mean([indices(1:end-1); indices(2:end)]));

indx = [indices(1) idx_for_integral indices(end)];

for k = 1:N
    num = indx(k + 1) - indx(k) + 1;
    S(k) = sum(S_(indx(k):indx(k+1)) .* dv)/(num .* dv); 
end

% plot(vel_axis_permitted,  abs(S)); hold on; plot(vel_axis,  abs(S_));
% 
% plot(S_(indices)); hold on; plot(S_(idx_for_integral));


% n = length(S); % Length of number of frequencies

idx = linspace(1, n, n);
idxq = linspace(min(idx), max(idx), N);

X = interp1(idx, X_full, idxq, 'linear');
Theta = interp1(idx, Theta_full, idxq, 'linear');
Noise = sum(S) ./ (N .* SNR); % Noise Power

P = -(S + Noise) .* log(X); % Power spectrum 

data_f = sqrt(P) .* exp(1j .* Theta); % frequency domain

% figure(1); hold on; plot(vel_axis, sqrt(abs(P))); grid on;

data = ifft(fftshift(sqrt(N) .* data_f));% complex time domain signal 
% data = ifft(fftshift(sqrt(P) .* exp(1j .* Theta)));% complex time domain signal 

% figure; plot(abs(data)); 

end