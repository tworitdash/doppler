% this one simulated Range-Doppler spectra map with and without noise


function [X_PS_nf,M0_truth,M1_truth,M2_truth] = weather_sig_simulator_beta(M0,M1,M2,v,dv)

L = numel(v); % number of sources

   
    ph = -pi+2*pi*rand(1,L);% uniformly distributed phase [-pi,pi]
    UD = rand(1,L);% uniformly distributed random variable [0,1]
    f = gauss_gen_2(v,M0,M1,M2);% v is in power
    
    X_PS_nf = ((-f.*log(UD))).*exp(1i*ph);% noise free Power Spectrum
      
    [M0_truth,M1_truth,M2_truth] = gauss_calc_2(abs(X_PS_nf),v,dv);% Moments computation with noise

end

% [1] D. S. Zrnic,Simulation of Weatherlike Doppler Spectra and Signals,J.Appl.Meteorol.14, no.4, 619 (June 1975)
% figure
% plot(10*log10(abs(X_PS_nf)),'r'),hold on
% plot(10*log10(abs(X_PS_n)),'b'),hold off