%% HD signal generator

close all;
clear;

% SNR_db = 30;
% SNR = 10^(SNR_db/10);
% lambda = 0.03;
% 
% x0 = 100;
% u = 2; 
% 
% NSweep = 5;
% n_rot = 10;
% NSec = 200;
% 
% Nt = n_rot * NSweep * NSec;
% 
% 
% dT = 1e-3; 
% x(1) = x0;
% 
% Z(1) = exp(1j * 4 * pi/lambda .* x(1));
% 
% for i = 2:Nt
%     x(i) = x(i - 1) + u * dT;
%     Z(i) = exp(1j * 4 * pi/lambda .* x(i)); 
% end
% 
% Noise = sum(abs(Z).^2)./(Nt .* SNR);
% sigma_n = sqrt(Noise);
% 
% Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2);
% 
% [~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);

%% Radar signal parameters with ground truth values 

r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 1.8*pi/180;
theta_0 = 0*pi/180;
dth = 2*pi/180;

% Scatterers positions

N = 2000;
r = r0-dr/2 + dr* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);
th = theta_0 - dth/2 + dth .* rand(1, N);

% positions of these scatteres in Cartesian domain

x0 = r .* cos(th) .* cos(ph);
y0 = r .* cos(th) .* sin(ph);
z0 = r .* sin(th); 

% r0 = sqrt(x0.^2 + y0.^2 + z0.^2);
r0 = sqrt(x0.^2 + y0.^2); % + z0.^2);

% figure; scatter3(x0, y0, z0, 50, r0); colormap('copper'); colorbar; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 

%% Signal model

D_min = 0.6;
D_max = 4;
D0 = 1;
N0 = 8e3;
u_mean = 1;
v_mean = 0;
u_sigma = 1;
v_sigma = 0;
SNR_db = 30;
n_rot = 1;
NSweep = 5;
NSec = 100;
Nt = n_rot * NSweep * NSec;
% Nt = 10;
dt = 1e-3;
lambda = 3e-2;

[Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma, sigma_n, Z] = ...
    Signal_modelvmac(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);

figure; histogram(real(Z)); figure; histogram(imag(Z));
figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));



%% HD Signal generator with mean and variance of U as parameters 

%% Available samples


Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep];
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dt;

figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

sigma_n_obs = sqrt(std(real(Z_avail_vec)));

figure; histogram(real(Z_model)); figure; histogram(real(Z_avail_vec));

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dt, lambda, SNR_db, 1, 1, 5);

[am, rm, as, rs, an, rn, itern] = MHradarOne([mu_obs sigma_obs sigma_n_obs], 10000000, (Z_avail_vec), t_avail, x0, y0, sigma_n);


burnin_num = round(length(am) .* 0.25);
am = am(burnin_num+1:end);
burnin_nums = round(length(as) .* 0.25);
as = as(burnin_nums+1:end);
burnin_numn = round(length(an) .* 0.25);
an = as(burnin_numn+1:end);


mu_re = mean(am); sigma_re = mean(as); sigma_n_re = mean(an);

figure; plot(rm, 'DisplayName', '\mu rejected'); hold on; plot(am, 'DisplayName', '\mu accepted'); legend;
figure; plot(rs, 'DisplayName', '\sigma rejected'); hold on; plot(as, 'DisplayName', '\sigma accepted'); legend;
figure; plot(rn, 'DisplayName', '\sigma_n accepted'); hold on; plot(an, 'DisplayName', '\sigma_n accepted'); legend;


figure; b = histogram(am); 

figure; hold on; plot(mu, '*', 'LineWidth', 5, 'DisplayName', '\mu GT'); 
hold on; plot(mu_re, '+', 'LineWidth', 5, 'DisplayName', '\mu MCMC mean'); 

[~, midx] = max(b.Values); mu_re_m = b.BinEdges(midx);
hold on; plot(mu_re_m, '^', 'LineWidth', 5, 'DisplayName', '\mu MCMC max'); legend;

figure; b = histogram(as);
figure; hold on; plot(sigma, '*', 'LineWidth', 5, 'DisplayName', '\sigma GT'); 
hold on; plot(sigma_re, '+', 'LineWidth', 5, 'DisplayName', '\sigma MCMC mean'); 

[~, midx] = max(b.Values); sigma_re_m = b.BinEdges(midx);
hold on; plot(sigma_re_m, '^', 'LineWidth', 5, 'DisplayName', '\sigma MCMC max');  legend;

figure; b = histogram(an);
figure; hold on; plot(sigma_n, '*', 'LineWidth', 5, 'DisplayName', '\sigma GT'); 
hold on; plot(sigma_n_re, '+', 'LineWidth', 5, 'DisplayName', '\sigma MCMC mean'); 

[~, midx] = max(b.Values); sigma_n_re_m = b.BinEdges(midx);
hold on; plot(sigma_n_re_m, '^', 'LineWidth', 5, 'DisplayName', '\sigma MCMC max');  legend;
%% Reconstruct
u_re = normrnd(mu_re, sigma_re, [1 N]);
t = eps:dt:(Nt - 1)*dt;
Z_re_(1) = sum(exp(1j .* 4 .* pi/lambda .* r0));
    
for l = 2:length(t)
    r = sqrt((x0 + u_re .* t(l)).^2 + y0.^2);
    Z_re_(l) = sum(exp(1j .* 4 .* pi/lambda .* r));
end

Z_re = Z_re_ + sigma_n .* (randn(1, Nt-1) + 1j .* randn(1, Nt-1))./sqrt(2);

[~, ~, mu_re_vr, sigma_re_vr, vel_axis_re, dv_re] = Spec(Z_re, Nt-1, dt, lambda, SNR_db, 1, 1, 4);