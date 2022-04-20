function ret = LLOne(x, data, t_avail, x0, sigma_n)

    lambda = 0.03;
    u = x;
    Z_model(1) = exp(1j .* 4 .* pi/lambda .* x0);
    
    for l = 2:length(t_avail)
        r = x0 + u .* t_avail(l);
        Z_model(l) = exp(1j .* 4 .* pi/lambda .* r);
    end
    
%     figure(1); plot(t_avail, real(Z_model)); hold on; plot(t_avail, imag(Z_model));
    
    ret_re = sum(-log(sigma_n .* sqrt(2 * pi)) - ((real(data)-real(Z_model).').^2)./(2 * sigma_n.^2));
    ret_im = sum(-log(sigma_n .* sqrt(2 * pi)) - ((imag(data)-imag(Z_model).').^2)./(2 * sigma_n.^2));
    
    ret = ret_re + 1j .* ret_im;
    

end