%% HD signal generator

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

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2);

[~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available samples


% Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
% Z_avail = Z_model_re(1:NSweep, :);
% Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

% figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));

Z_model_re = reshape(Z, [NSweep * NSec n_rot]);
Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]);

Noise_avail = sum(abs(Z_avail_vec_).^2)./(length(Z_avail_vec_) .* SNR);
sigma_n_avail = sqrt(Noise_avail);

Z_avail_vec = Z_avail_vec_ + sigma_n_avail .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep];
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT;


% figure; plot(t_avail, real(Z_avail_vec_)); hold on; plot(t_avail, real(Z_m), '*');

% figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);


%% MCMC parameters in a structure var
sigma_n_0 = 0.1;

% E.n = 2;
% E.E0 = [mu_obs sigma_n_0];
% E.sig = [sigma_obs-1.5 0.05];
% E.H = [mu_obs + (sigma_obs-1.5) 0.5];
% E.L = [mu_obs - (sigma_obs - 1.5) 0.01];
% [accepted, rejected, itern] = MHu(E, 10000, Z_avail_vec, t_avail, x0);
E.n = 1;

E.E0 = [mu_obs];
E.sig = [50000];
E.H = [mu_obs + (sigma_obs-1.5)];
E.L = [mu_obs - (sigma_obs - 1.5)];
[accepted, rejected, itern, E_new] = MHu(E, 1000, Z_avail_vec, t_avail, x0, sigma_n);


for i = 1:E.n

    
   figure(1000+i);plot(rejected(:, i)); hold on; plot(accepted(:, i));
   
%    Mest(i) = mean(accepted(:, i));
   burnin = round(0.25 * length(accepted(:, i)));
   figure(2000+i); histogram(accepted(burnin+1:end, i), 100);
   burninrej = round(0.25 * length(rejected(:, i)));
   figure(3000+i); histogram(rejected(burninrej+1:end, i), 100);
   mu_re = mean(accepted(burnin+1:end, i)); 
   rej_re = mean(rejected(burnin+1:end, i));
end


% sigma_n_re = Mest(2);





