function ll = like(u, z, t, sigma_n, mu_0, sigma_2u)
    lambda = 0.03; 
    x = exp(1j .* u * t*4*pi/lambda);
    n = z - x;
    n_re_std = real(n) * sqrt(2) /sigma_n;
    n_im_std = imag(n) * sqrt(2) /sigma_n;
    
    ll1 = -1/2*log(2*pi) - 1/2 * n_re_std.^2;
    ll2 = -1/2*log(2*pi) - 1/2 * n_im_std.^2;
    
    log_pu = -1/2*log(2*pi*sigma_2u)-1/2*(u-mu_0)^2/sigma_2u;
    
    ll = sum(ll1) + sum(ll2)  + log_pu;
end