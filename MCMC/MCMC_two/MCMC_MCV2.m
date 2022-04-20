%% HD signal generator

close all;
clear;

SNR_db = 30;
SNR = 10^(SNR_db/10);
lambda = 0.03;

x0 = 100;
u = 4; 

NSweep = 5;
n_rot = 10;
NSec = 200;

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

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2);

[~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available samples


% Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
% Z_avail = Z_model_re(1:NSweep, :);
% Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

% figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));

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

figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);

MC = 16; 

sigma_n_obs = std(real(Z_avail_vec));

%% MCMC parameters in a structure var

E.n = 2;
E.E0 = [mu_obs sigma_n_obs];
E.sig = [sigma_obs 2];
E.H = [mu_obs + sigma_obs sigma_n_obs];
E.L = [mu_obs - sigma_obs 0];

%% Running the MCMC chain 

for l = 1:MC

        [a, r, itern] = MH(E, 10000, (Z_avail_vec), t_avail, x0);

    
        burnin_num = round(size(a, 1) .* 0.25);
        
        a = a(burnin_num+1:end, :);
        re(l, :) = mean(a, 1);
        
        figure(1000); plot(a);
        
  
end

 mu_re = mean(re(:, 1), 1);
 sigma_n_re = mean(re(:, 2), 1);


% figure; plot(r); hold on; plot(a);

% figure; b = histogram(a);
% [~, midx] = max(b.Values); mu_re_m = b.BinEdges(midx);


%% Reconstruct
t = eps:dT:(Nt - 1)*dT;
Z_re_(1) = exp(1j .* 4 .* pi/lambda .* x0);
    
for l = 2:length(t)
    r = x0 + mu_re .* t(l);
    Z_re_(l) = exp(1j .* 4 .* pi/lambda .* r);
end

Z_re = Z_re_ + sigma_n_re .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2);

[~, ~, mu_re_vr, sigma_re_vr, vel_axis_re, dv_re] = Spec(Z_re, Nt, dT, lambda, SNR_db, 1, 1, 4);