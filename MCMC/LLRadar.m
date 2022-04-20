function ret = LLRadar(x, data, t_avail, x0, y0, N, sigma_n)
    lambda = 0.03;
    U = normrnd(x(1), x(2), [1 N]);
    
    
    for l = 1:length(t_avail)
        r = sqrt((x0 + U .* t_avail(l)).^2 + y0.^2);
        Z_model(l) = sum(exp(1j .* 4 .* pi/lambda .* r));
    end
    
    ret_re = sum(-log(sigma_n .* sqrt(2 * pi)) - ((real(data)-real(Z_model)).^2)./(2 * sigma_n.^2));
    ret_im = sum(-log(sigma_n .* sqrt(2 * pi)) - ((imag(data)-imag(Z_model)).^2)./(2 * sigma_n.^2));
    
    ret = ret_re + 1j .* ret_im;
end