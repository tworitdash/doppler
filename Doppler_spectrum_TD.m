%% Function to calculate the time domain signal for radar when the maximum and minimum Doppler velocities are given as inputs

function [data, sig_f] = Doppler_spectrum_TD(vel_axis, mu, sigma, Nifft, SNR, lambda, X, Theta)

% freq_axis = 2 .* vel_axis ./ lambda; % Doppler frequency axis

% S = zeros(size(freq_axis)); % Signal dependent on frequency

% mu = (v_min + v_max)./2; % Mean Doppler velocity needed
% sigma = (v_max - v_min)./10;

% [~, vel_indices_required] = find(vel_axis>v_min&vel_axis<v_max); % Find the indices and the vel spectrum of interest

% S(vel_indices_required) = ones(1, length(vel_indices_required)); % normpdf(vel_axis(vel_indices_required), mu, sigma);  % Signal Spectrum (frequency dependent)

S = 1/sqrt(2*pi*sigma^2) * exp(-(vel_axis - mu).^2/(2*sigma^2));  

n = length(S); % Length of number of frequencies

N = sum(S) ./ (n .* SNR); % Noise Power

P = -(S + N) .* log(X); % Power spectrum 
sig_f = sqrt(P);

% figure(1); hold on; plot(vel_axis, sqrt(abs(P))); grid on;

data = ifft(fftshift(sqrt(Nifft) .* sqrt(P) .* exp(1j .* Theta)));% complex time domain signal 
% data = ifft(fftshift(sqrt(P) .* exp(1j .* Theta)));% complex time domain signal 

end