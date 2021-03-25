function [data] = data_simulator_BW(time_axis, beta_i, beta_wind, SNR, lambda, vel)
    data = zeros(1, length(time_axis));
    noise = zeros(1, length(time_axis));
    for i = 1:length(time_axis)
        beta = (beta_wind - beta_i);     % Difference between wind direction and radar beam 
        vr = vel .* cos(beta);
        noise(i) = 1/(SNR) .* (randn(1) + 1j .* randn(1));
        
%         data(i) = exp(1j .* 2.*pi/lambda .* 2 .* (vel) .* time_axis(i)) + noise(i); %1/sqrt(SNR) .* randn(1);
        data(i) = exp(1j .* 2.*pi/lambda .* 2 .* (vr) .* time_axis(i)) + noise(i); %1/sqrt
    end
%     data_doppler = 1./sqrt(length(data)) .* fftshift(fft(data));
%     Theta = rand(1, length(data)) .* 2 * pi; 
% %     Theta = 0;
%     data_ = ifft(fftshift(sqrt(length(data)) .* (data_doppler) .* exp(1j .* Theta)));
%     data = abs(data_) .* exp(1j .* unwrap(angle(data_)) .* cos(beta));
end