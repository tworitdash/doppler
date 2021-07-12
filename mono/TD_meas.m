function [s] = TD_meas(mu_0, lambda, beta_wind_0, phi_0, Omega, t)

    alpha_wind = pi/3/t(end);
    
    mu_der = mu_0./t(end);
    
    C = beta_wind_0 - phi_0;
    D = alpha_wind - Omega;
    
    ph_int = (D .* (mu_0 + mu_der .* t) .* sin(C + D .* t) + mu_der .* cos(C + D .* t)) ./ D.^2;
    
    s = exp(1j .* 4 .* pi ./ lambda .* ph_int);
  
end