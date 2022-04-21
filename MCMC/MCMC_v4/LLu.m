function ret = LLu(u, data, t_avail, r0, sigma_n, pu)

    lambda = 0.03;
    
    r(1) = r0; % initial position of the scatterer
    Z_model(1) = exp(1j .* 4 .* pi/lambda .* r(1)); % First sample echo model 
    
    for l = 2:length(t_avail)
            r(l) = r(l - 1) + u .* (t_avail(l) - t_avail(l - 1)); 
            % This for loop calculates the model
            % samples from the available t steps $t_avail$
          Z_model(l) = exp(1j .* 4 .* pi/lambda .* r(l));
    end
    
    Nt = length(data);
   
    K = eye(Nt, Nt) .* sigma_n.^2; % Noise diagnoal matrix with size length(t_avail) x length(t_avail)
   

%% Real and imaginary part likelihood 

    ret_re = log(pu) -Nt/2 .* log(2*pi) - 1/2 * (real(data) - real(Z_model).').' * inv(K) * (real(data) - real(Z_model).');
    ret_im = log(pu) -Nt/2 .* log(2*pi) - 1/2 * (imag(data) - imag(Z_model).').' * inv(K) * (imag(data) - imag(Z_model).');
    
    ret = ret_re + ret_im; % Total likelihood (multiplication of real and imaginary likelihoods in linear scale)
    

end