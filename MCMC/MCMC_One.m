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


Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep];
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT;

figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);

[a, r, itern] = MHradarOne([mu_obs sigma_obs], 10000, (Z_avail_vec), t_avail, x0, sigma_n);


burnin_num = round(length(a) .* 0.25);
a = a(burnin_num+1:end);
mu_re = mean(a);

figure; plot(r); hold on; plot(a);

figure; b = histogram(a);
[~, midx] = max(b.Values); mu_re_m = b.BinEdges(midx);


%% Reconstruct
t = eps:dT:(Nt - 1)*dT;
Z_re_(1) = exp(1j .* 4 .* pi/lambda .* x0);
    
for l = 2:length(t)
    r = x0 + mu_re .* t(l);
    Z_re_(l) = exp(1j .* 4 .* pi/lambda .* r);
end

Z_re = Z_re_ + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2);

[~, ~, mu_re_vr, sigma_re_vr, vel_axis_re, dv_re] = Spec(Z_re, Nt, dT, lambda, SNR_db, 1, 1, 4);