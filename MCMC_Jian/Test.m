close all;
clear;

SNR_db = 80;
SNR = 10^(SNR_db/10);
lambda = 0.03;

x0 = 100;
u = 4; 

NSweep = 5;
n_rot = 100;
NSec = 20;

Nt = n_rot * NSweep * NSec;


dT = 1e-3; 
x(1) = x0;

Z(1) = exp(1j * 4 * pi/lambda .* x(1));

for i = 2:Nt
    x(i) = x(i - 1) + u * dT;
    Z(i) = exp(1j * 4 * pi/lambda .* x(i)); 
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);
sigma_n = sqrt(Noise);

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt));
   
Z_model_re = reshape(Z, [NSweep * NSec n_rot]);
Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

Noise_avail = sum(abs(Z_avail_vec).^2)./(length(Z_avail_vec) .* SNR);
sigma_n_avail = sqrt(Noise_avail);

Z_avail_vec = Z_avail_vec + sigma_n_avail .* (randn(1, length(Z_avail_vec)).'+ 1j .* randn(1, length(Z_avail_vec)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep];
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; 

%% 
%     V = linspace(0, 6, 1000);
%     
%     for i = 1:1000
        
    u1 = 1;
    Z_u1(1) = exp(1j .* 4 .* pi/lambda .* x0);
    
    for l = 2:length(t_avail)
        r = x0 + u1 .* t_avail(l);
        Z_u1(l) = exp(1j .* 4 .* pi/lambda .* r);
    end
    
    N_avail = length(Z_avail_vec);
    
    K = eye(N_avail, N_avail) .* sigma_n.^2;
    
%     detK = x(2).^(2*Nt);
%     detK = sigma_n.^(2*Nt);
    
%     ret_re = -Nt/2 .* log(2*pi) - 1/2 .* log(detK) - 1/2 * (real(data) - real(Z_model).').' * inv(K) * (real(data) - real(Z_model).');
%     ret_im = -Nt/2 .* log(2*pi) - 1/2 .* log(detK) - 1/2 * (imag(data) - imag(Z_model).').' * inv(K) * (imag(data) - imag(Z_model).');
     
    ret_re = - 1/2 * (real(Z_avail_vec) - real(Z_u1).').' * inv(K) * (real(Z_avail_vec) - real(Z_u1).');
    ret_im = - 1/2 * (imag(Z_avail_vec) - imag(Z_u1).').' * inv(K) * (imag(Z_avail_vec) - imag(Z_u1).');

    
%     ret_re = sum(-log(sigma_n .* sqrt(2 * pi)) - ((real(data)-real(Z_model).').^2)./(2 * sigma_n.^2));
%     ret_im = sum(-log(sigma_n .* sqrt(2 * pi)) - ((imag(data)-imag(Z_model).').^2)./(2 * sigma_n.^2));
    
    ret = ret_re + 1j .* ret_im;
    
%     an(i) = unwrap(angle(ret));
%   
%    
%     end
    
%     figure; plot(V, an);