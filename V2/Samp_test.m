close all;
clear;
SNR = 10^(30/10);
m0 = 1;
mu = 5;
sigma = 0.2;
n = 2^10;
v_amb = 7.5;
PRT = 1e-3;
lambda = 0.03;

N  = 2.^(linspace(1, 5, 5));

for i = 1:length(N)
    [dataB, data_fB, ~, ~, X_full, Theta_full, vel_axisB, dvB] = DS_simulatorV2B(SNR, m0, mu, sigma, n, v_amb, N(i));
    data_fB2 = 1/sqrt(N(i) + 1) .* fftshift(fft(dataB, N(i) + 1));
  
    
    
    PTB = sum(abs(data_fB2).^2 .* dvB);
    numB = sum(vel_axisB .* abs(data_fB2).^2 .* dvB);
    vB(i) = numB./PTB;
%     
    
    [dataA, data_fA, ~, ~, ~, ~, vel_axisA, dvA] = DS_simulatorV2A(SNR, m0, mu, sigma, n, v_amb, N(i), X_full(1:N(i)), Theta_full(1:N(i)));
    data_fA2 = 1/sqrt(N(i)) .* fftshift(fft(dataA, N(i)));
    
    PTA = sum(abs(data_fA2).^2 .* dvA);
    numA = sum(vel_axisA .* abs(data_fA2).^2 .* dvA);
    vA(i) = numA./PTA;
    
    
    data_fA2_fft = 1/sqrt(N(i)) .* (fft(dataA, N(i)));
    
    [~, indx] = max(abs(data_fA2_fft).^2);
    
    
    PT_fft = sum(abs(data_fA2_fft).^2); % Total power of the Doppler Spectrum
             
    axis_mean = indx-N(i)/2:1:indx+N(i)/2-1;
                
    f_mean_l(i) = (1/(N(i) .* PRT)) .* indx + (1/(N(i) .* PRT .* PT_fft))...
                    .* sum((axis_mean - indx) .* abs(data_fA2_fft).^2);
    vA_fft(i) = f_mean_l(i) .* lambda./2;
                
%     if vA_fft(i) > v_amb
%         vA_fft(i) = vA_fft(i) - 2.*v_amb;
%     end
%     
    
   
    figure(101); hold on; plot(vel_axisA, (abs(data_fA))); grid on; 
    
    hold on; plot(vel_axisA, (abs(data_fA2)), 'o')
%     
%   
%     
    figure(102); hold on; plot(vel_axisB, (abs(data_fB))); grid on;
    hold on; plot(vel_axisB, (abs(data_fB2)), 'o')
end

figure; plot(N, vA); hold on; plot(N, vB, '*-'); % hold on; plot(N, vA_fft, '*'); grid on; 

%%
mu = 5;
v_amb = 7.5;
lambda = 0.03;
PRT = 1e-3;
% N = 2.^(linspace(2, 5, 4));
N = 8;


for i = 1:length(N(1))
    sig = exp(1j .* 2 .* pi .* 2 .* mu / lambda .* (1:N(i)) .* PRT);
    dv = 2*v_amb/N(i);
    sig_f = 1/sqrt(N(i)) .* fftshift(fft(dataB, N(i)));
    vel_axis = linspace(-v_amb, v_amb-dv, N(i));
    sig_f2 = 1/sqrt(N(i) + 1) .* fftshift(fft(dataB, N(i) + 1));
    vel_axis2 = linspace(-v_amb, v_amb, N(i)+1);
    figure;
    plot(vel_axis, abs(sig_f).^2);
    hold on;
    plot(vel_axis2, abs(sig_f2).^2); grid on;
end
