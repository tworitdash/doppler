function ret = LL(x, data, t_avail, x0, y0)

    lambda = 0.03;
    N = 200;
    u = normrnd(x(1), x(2), [1 N]);
    Z_model(1) = sum(exp(1j .* 4 .* pi/lambda .* sqrt(x0.^2 + y0.^2)));
    
    for l = 2:length(t_avail)
        r = sqrt((x0 + u .* t_avail(l)).^2 + y0.^2);
        Z_model(l) = sum(exp(1j .* 4 .* pi/lambda .* r));
    end
    
%     figure(1); plot(t_avail, real(Z_model)); hold on; plot(t_avail, imag(Z_model));
    
    ret_re = sum(-log(x(3) .* sqrt(2 * pi)) - ((real(data)-real(Z_model).').^2)./(2 * x(3).^2));
    ret_im = sum(-log(x(3) .* sqrt(2 * pi)) - ((imag(data)-imag(Z_model).').^2)./(2 * x(3).^2));
    
    ret = ret_re + 1j .* ret_im;
    

end