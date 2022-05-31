%% 
clear; 
% close all;

%% Scatterer information initial 

N = 2000; % Number of scatterers

r0 = 100;
phi_deg = 0;
phi0 = phi_deg * pi/180;
dr = 50;
dphi_deg = 2;
dphi = dphi_deg * pi/180;

raxis(1, :) = r0 + (dr) .* rand(1, N);
phiaxis(1, :) = phi0 + (dphi) .* rand(1, N);

r0_ = linspace(r0, r0 + dr, N);
phi0_ = linspace(phi0, phi0 + dphi, N);

[r0_, phi0_] = meshgrid(r0_, phi0_);


x0 = r0_ .* cos(phi0_);
y0 = r0_ .* sin(phi0_);

x(1, :) = raxis(1, :) .* cos(phiaxis(1, :));
y(1, :) = raxis(1, :) .* sin(phiaxis(1, :));
r(1, :) = sqrt(x(1, :).^2 + y(1, :).^2);

%% Velocity profiles
MU_U = 3;
W_U = 0.1;

MV_V = 0;
W_V = 0;

U = normrnd(MU_U, W_U, [1 N]); 
V = normrnd(MV_V, W_V, [1 N]);

%% Time dimension details

lambda = 3e-2;
dT = 1e-3;
NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
n_rot = 1;                % Number of rotations of radar
NSec = 200;                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)

t = linspace(eps, (Nt - 1).*dT, Nt);

%% Plot the phsyical model 

% figure; pcolor(x0, y0, ones(N, N)); % shading flat; 
% hold on;
% scatter(x(1, :),y(1, :),40,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0.6350 .0780 .1840],...
%               'LineWidth',1.5)

for k = 2:Nt
    x(k, :) = x(k-1, :) + U .* dT;
    y(k, :) = y(k-1, :) + V .* dT;
%     pause(0.1);
%     pcolor(x0, y0, ones(N, N)); shading flat; 
%     hold on;
%     scatter(x(k, :),y(k, :),40,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0.6350 .0780 .1840],...
%               'LineWidth',1.5)
end

%% Radar echoes from the model


SNR_db = 30;
SNR = 10^(SNR_db/10);
v_amb = lambda/(2 * dT);

s(1) = sum(exp(1j .* 4 * pi / lambda * r(1, :)));

for k = 2:Nt
    r(k, :) = sqrt(x(k, :).^2 + y(k, :).^2);
    s(k) = sum(exp(1j .* 4 * pi / lambda * r(k, :)));
    
    vr(k-1, :) = (r(k, :) - r(k-1, :))/dT;
end

% figure(1); hold on; plot(real(s)); figure(2); hold on; plot(imag(s));

Noise = sum(abs(s).^2)./(Nt .* SNR);
n = sqrt(Noise)/sqrt(2) .* (randn(1, Nt) + 1j .* randn(1, Nt));
z = s + n;

[~, ~, vrmu, vrsigma, vel_axis, dv] = Spec(z, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available measurement model 



% Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
% Z_avail = Z_model_re(1:NSweep, :);
% Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

% figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));

Z_model_re = reshape(z, [NSweep * NSec n_rot]); 
Z_avail = Z_model_re(1:NSweep, :);              % Takes only N_Sweep number of samples from all rotations             
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]); % vectorize all samples available

% Noise_avail = sum(abs(Z_avail_vec_).^2)./(length(Z_avail_vec_) .* SNR);
% sigma_n_avail = sqrt(Noise_avail);

Z_avail_vec = Z_avail_vec_; %  + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t_(:, k) = (k - 1) * NSec * NSweep + [0:NSweep-1]; % This for loop calculates the time instances of the available samples
end

t_avail = reshape(t_, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances

N_avail = length(t_avail);

% figure; plot(t_avail, real(Z_avail_vec_)); hold on; plot(t_avail, real(Z_m), '*');

% figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

% This functions calculates the Doppler spectrum of the available samples
% in one rotation only N_Sweeps coherent samples. mu_obs and sigma_obs are
% the mean Doppler velocity and Doppler spectrum width with only N_Sweeps
% samples. Mu_obs is taken as the starting point for the MCMC algorithm.

[~, ~, vrmu_obs, vrsigma_obs, vel_axis_obs, dv_obs] =  ...
    Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);

% [~, ~, vrmu_obs, vrsigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec.', N_avail, dT, lambda, SNR_db, 1, 1, 5);

%% MCMC algorithm

MCMC.n = 2;
MCMC.E0 = [vrmu_obs vrsigma_obs];
MCMC.H = [vrmu_obs+2, 3];
MCMC.L = [0, 0];

Num_Simu = 100000;

% Known quantities 

[accepted, rejected, itern, E_new, sample_MC_seq, sample_all] ...
= MHvr2D(MCMC, Num_Simu, Z_avail_vec, t_avail, N, x(1, :), y(1, :), sqrt(Noise));