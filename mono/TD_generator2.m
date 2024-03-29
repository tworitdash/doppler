function [s, SNR] = TD_generator2(mu, lambda, beta_wind, phi_0, Omega, t, n_sig)
% 
%     ph_ = (4 * pi/lambda * mu .* t);
%     s = (exp(1j .* ph_ .* (sin(eps + beta_wind - Omega .* t - phi_0)./(eps - Omega .* t))));
    
    PRT = t(2) - t(1);
    
%     figure; plot(mu);
    
    vel = mu .* cos(beta_wind - Omega .* t - phi_0);
    
    for i = 1:length(t)
        dis(i) = sum(vel(1:i) .* PRT);
    end
    
    ph_ = 4 * pi / lambda .* dis;
    s = exp(1j .* ph_) + n_sig .* randn(size(t));
    
    SNR = sum(abs(s).^2)/(length(t) .* n_sig^2);

end