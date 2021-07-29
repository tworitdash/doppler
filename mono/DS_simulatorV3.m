function [data, data_f, data_f_Sig, X, Theta, P_wNoise, P] = DS_simulatorV3(SNR, m0, mu, sigma, N, v_amb, p)

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


% if sigma < 0.02
%     [~, idx1] = (min(abs(vel_axis - mu)));
%     S_ = dirac(vel_axis - vel_axis(idx1)); 
%     idx = S_ == Inf;
%     S_(idx) = 1;
% else
for i = 1:N
if vel_axis(i) < mu
    S_(i) = (2/(1 + p)) .* m0./sqrt(2*pi*(sigma).^2) .* exp(-(vel_axis(i) - mu).^2./(2*(sigma).^2));
else 
    S_(i) = (2 .* p./ (1 + p)) .* m0./sqrt(2*pi*(sigma .* p).^2) .* exp(-(vel_axis(i) - mu).^2./(2*(sigma .* p).^2));
end
end

% figure; plot(abs(S_))
% end
%% Doviak

% Noise = sum(S_) ./ (N .* SNR); % Noise Power
% 
% P_wNoise = -(S_ + Noise) .* log(X); % Power spectrum 
% 
% 
% % P = -S_;
% data_f = sqrt(P_wNoise) .* exp(1j .* Theta);
% data = ifft(ifftshift(sqrt(N) .* data_f));
% %% Only signal and no noise
% P = -(S_); % Power spectrum 
% 
% data_f_Sig = sqrt(P) .* exp(1j .* Theta);

%% Sirmans

PN = 1;

K = (PN .* SNR)./(sum(S_));

P_wNoise = -log(X) .* (K .* S_ + PN./N);

data_f = sqrt(P_wNoise) .* exp(1j .* Theta);
data = ifft(ifftshift(sqrt(N) .* data_f));
%% Only signal and no noise
P = -(S_); % Power spectrum 

data_f_Sig = sqrt(P) .* exp(1j .* Theta);


end