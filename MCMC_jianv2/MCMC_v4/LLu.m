function ret = LLu(x, data, t_avail, x0, sigma_n, px)

    lambda = 0.03;
    c = 3e8;
    u = x(1); % u is chosen 
    
    r(1) = x0; % initial position of the scatterer
    Z_model(1) = exp(1j .* 4 .* pi/lambda .* r(1)); % First sample echo model 
    
    for l = 2:length(t_avail)
            r(l) = r(l - 1) + u .* (t_avail(l) - t_avail(l - 1)); 
            % This for loop calculates the echo...
            % samples from the available time steps $t_avail$
          Z_model(l) = exp(1j .* 4 .* pi/lambda .* r(l));
    end
    Nt = length(data);
 
    K = eye(2*Nt, 2*Nt) .* sigma_n.^2; % Noise diagnoal matrix with size length(t_avail) x length(t_avail)
    data_ri = [real(data); imag(data)];
    Z_model_ri = [real(Z_model), imag(Z_model)].';

    ret = log(px) -Nt * log(2*pi) - Nt*2*log(sigma_n) - 1/(2*sigma_n^2) * (data_ri - Z_model_ri).'  * (data_ri - Z_model_ri) ;

end