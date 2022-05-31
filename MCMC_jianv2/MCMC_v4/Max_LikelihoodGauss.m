function ret = Max_LikelihoodGauss(x, Nscatters, data, t_avail, x0, sigma_n, px)
%     Nscatters = x(Nscatters + 1);
    r1(1, :) = x0 .* rand(1, Nscatters);
%     r1(1, :) = x0;

    lambda = 0.03;
    c = 3e8;
%     a  = x(3);
    a = 1;
    u1 = x(1) - x(2)/2 + x(2) .* rand(1, Nscatters); % u is chosen 
%     u2 = x(2);
    
%     r1(1) = x0; % initial position of the scatterer
%     r2(1) = x0;
    Z_model(1) = sum( exp(1j .* 4 .* pi/lambda .* r1(1, :)) ) ; % First sample echo model 
    
    for l = 2:length(t_avail)
           r1(l, :) = r1(l - 1, :) + u1 .* (t_avail(l) - t_avail(l - 1));
%            r2(l) = r2(l - 1) + u2 .* (t_avail(l) - t_avail(l - 1)); 
            % This for loop calculates the echo...
            % samples from the available time steps $t_avail$
           Z_model(l) = sum( exp(1j .* 4 .* pi/lambda .* r1(l, :)) );
    end
    Nt = length(data);
 
    K = eye(2*Nt, 2*Nt) .* sigma_n.^2; % Noise diagnoal matrix with size length(t_avail) x length(t_avail)
    data_ri = [real(data); imag(data)];
    Z_model_ri = [real(Z_model), imag(Z_model)].';

    ret = abs(log(px) - Nt * log(2*pi) - Nt*2*log(sigma_n) - 1/(2*sigma_n^2) * (data_ri - Z_model_ri).'  * (data_ri - Z_model_ri)) ;

end