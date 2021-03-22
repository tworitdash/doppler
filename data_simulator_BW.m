function [data] = data_simulator_BW(time_axis, beta_i, beta_wind, SNR, lambda, vel)
    data = zeros(1, length(time_axis));
    noise = zeros(1, length(time_axis));
    for i = 1:length(time_axis)
        beta = (beta_wind - beta_i);     % Difference between wind direction and radar beam 
        vr = vel .* cos(beta);
        noise(i) = 1/(SNR) .* (randn(1) + 1j .* randn(1));
        data(i) = exp(1j .* 2.*pi/lambda .* 2 .* (vr) .* time_axis(i)) + noise(i); %1/sqrt(SNR) .* randn(1);
    end
end