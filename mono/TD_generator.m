function [s] = TD_generator(mu, lambda, beta_wind, phi_0, Omega, t)

    ph_ = (4 * pi/lambda * mu .* t);
    s = (exp(1j .* ph_ .* (sin(eps + beta_wind - Omega .* t - phi_0)./(eps - Omega .* t))));

end