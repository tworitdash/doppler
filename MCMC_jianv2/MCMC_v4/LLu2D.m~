function ret = LLu2D(U, data, t_avail, x0, y0, sigma_n, px)

    lambda = 0.03;
    c = 3e8;
    u = U(1); % u is chosen 
    v = U(2); % v is chosen
    B = U(3);
    x(1) = x0;
    y(1) = y0;
    
    r(1) = sqrt(x0.^2 + y0.^2); % initial position of the scatterer
    Z_model(1) = B .* exp(1j .* 4 .* pi/lambda .* r(1)); % First sample echo model 
    
    for l = 2:length(t_avail)
            x(l) = x(l - 1) + u .* (t_avail(l) - t_avail(l - 1));
            y(l) = y(l - 1) + v .* (t_avail(l) - t_avail(l - 1));
            r(l) = sqrt((x(l)).^2 + (y(l)).^2); 
            % This for loop calculates the echo...
            % samples from the available time steps $t_avail$
            Z_model(l) = B .* exp(1j .* 4 .* pi/lambda .* r(l));
    end
    Nt = length(data);
 
    K = eye(2*Nt, 2*Nt) .* sigma_n.^2; % Noise diagnoal matrix with size length(t_avail) x length(t_avail)
    data_ri = [real(data); imag(data)];
    Z_model_ri = [real(Z_model), imag(Z_model)].';

    ret = log(px) - Nt * log(2*pi) - Nt*2*log(sigma_n) - 1/(2*sigma_n^2) * (data_ri - Z_model_ri).'  * (data_ri - Z_model_ri) ;

end