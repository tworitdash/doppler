function ret = LLu(x, data, t_avail, x0, sigma_n, px)

    lambda = 0.03;
    u = x(1); % u is chosen 
    
    r(1) = x0; % initial position of the scatterer
    Z_model(1) = exp(1j .* 4 .* pi/lambda .* r(1)); % First sample echo model 
    
    for l = 2:length(t_avail)
            r(l) = r(l - 1) + u .* (t_avail(l) - t_avail(l - 1)); 
            % This for loop calculates the echo...
            % samples from the available time steps $t_avail$
          Z_model(l) = exp(1j .* 4 .* pi/lambda .* r(l));
    end
    
%     for l = 2:length(t_avail)
%         r = x0 + u .* t_avail(l);
%         Z_model(l) = exp(1j .* 4 .* pi/lambda .* r);
%     end
    
%     figure(1); plot(t_avail, real(Z_model)); hold on; plot(t_avail, imag(Z_model));
    Nt = length(data);
    
%     K = eye(Nt, Nt) .* x(2).^2;
    K = eye(Nt, Nt) .* sigma_n.^2; % Noise diagnoal matrix with size length(t_avail) x length(t_avail)
    
%     detK = x(2).^(2*Nt);
%     detK = sigma_n.^(2*Nt);
    
%     ret_re = -Nt/2 .* log(2*pi) - 1/2 .* log(detK) - 1/2 * (real(data) - real(Z_model).').' * inv(K) * (real(data) - real(Z_model).');
%     ret_im = -Nt/2 .* log(2*pi) - 1/2 .* log(detK) - 1/2 * (imag(data) - imag(Z_model).').' * inv(K) * (imag(data) - imag(Z_model).');
     

%% Real and imaginary part likelihood 

    ret_re = log(px) -Nt/2 .* log(2*pi) - 1/2 * (real(data) - real(Z_model).').' * inv(K) * (real(data) - real(Z_model).');
    ret_im = log(px) -Nt/2 .* log(2*pi) - 1/2 * (imag(data) - imag(Z_model).').' * inv(K) * (imag(data) - imag(Z_model).');

    
%     ret_re = sum(-log(sigma_n .* sqrt(2 * pi)) - ((real(data)-real(Z_model).').^2)./(2 * sigma_n.^2));
%     ret_im = sum(-log(sigma_n .* sqrt(2 * pi)) - ((imag(data)-imag(Z_model).').^2)./(2 * sigma_n.^2));
    
    ret = ret_re + 1j .* ret_im;
    

end