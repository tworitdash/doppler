%% HD signal generator

close all;
clear;

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
n_rot = 100;
NSweep = 5;
NSec = 100;
Nt = n_rot * NSweep * NSec;
% Nt = 10;
dt = 1e-3;
lambda = 3e-2;

[Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma, sigma_n, Z] = ...
    Signal_modelvmac(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db);
% 
% figure; histogram(real(Z)); figure; histogram(imag(Z));
figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));

mu_gtr = mean(real(Z_model));
sigma_gtr = std(real(Z_model));

mu_gti = mean(imag(Z_model));
sigma_gti = std(imag(Z_model));

GT_param = [mu_gtr sigma_gtr mu_gti sigma_gti];

E.gt = GT_param;

%% HD Signal generator with mean and variance of U as parameters 

%% Available samples


Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
Z_avail = Z_model_re(1:NSweep, :);
Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep];
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dt;

% figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))
figure; histogram(real(Z_avail_vec)); figure; histogram(imag(Z_avail_vec));

mu_sig_obsr = mean(real(Z_avail_vec));
std_obsr = std(real(Z_avail_vec));

mu_sig_obsi = mean(imag(Z_avail_vec));
std_obsi = std(imag(Z_avail_vec));

%% MCMC on the signal properties 


E.n = 4;
E.E0 = [mu_sig_obsr std_obsr mu_sig_obsi std_obsi];
E.sig = [std_obsr N/10 std_obsi N/10];
E.H = [mu_sig_obsr + std_obsr N mu_sig_obsi + std_obsi N] ;
E.L = [mu_sig_obsr - std_obsr 0 mu_sig_obsi - std_obsi 0 ];

ii = ["\mu Real part" "\sigma Real part" "\mu Imaginary part" "\sigma Imaginary part"];

%% 

[a, r, itern] = MH(E, 100, (Z_avail_vec), t_avail, x0);

burnin_num = round(size(a, 1) .* 0.25);
        
a = a(burnin_num+1:end, :);
% re = mean(a, 1);


for i = 1:E.n
    figure(i+100); plot(r(:, i)); hold on; plot(a(:, i),  'DisplayName', ii(i)); 
    figure(i+200); histogram(a(:, i), 'DisplayName', ii(i));
end

%% Reconstruct 

mean_val = mean(a, 1);

Z_re = normrnd(mean_val(1), mean_val(2), [1 Nt]) + 1j .* normrnd(mean_val(3), mean_val(4), [1 Nt]);

figure(1); hold on; histogram(real(Z_re)); figure(2); hold on; histogram(imag(Z_re));

% Spec(Z_re, Nt, dt, lambda, SNR_db, 1, 1, 7);

%% interactive check 

% 
% fig1 = figure('Position',[0 0 800 800]);
% fig2 = figure('Position',[850 0 800 800]);
% 
% figure(fig2); histogram(real(Z_model));
% 
% for k = 1:Nt
% %     figure(fig1); hold on; plot(a(k, 1));
%     figure(fig2); histogram(real(Z_re(1:k))); hold off;
% end
