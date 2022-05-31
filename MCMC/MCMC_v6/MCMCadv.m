clear; 
close all;

lambda = 0.03; 
dt = 1e-3;

u = 4;
SNR_db = 30;
SNR = 10^(SNR_db/10);

n_rot = 10;
K = linspace(0, n_rot-1, n_rot);

M = 5;
Nsec = 200;
N = Nsec * M;

Nt = N * n_rot;

t = linspace(0, (Nt -1)*dt, Nt);
t_avail = [];
for k = 0:n_rot-1
    for i = 0:M-1
        t_avail = [t_avail N * (k) * dt + (i) * dt];
    end
end
%% Target info

r0 = 0;

%% HD signal

xg = exp(1j * (r0 + u * t) * 4 * pi/lambda);

% Noiseg = sum(abs(xg).^2)./(Nt .* SNR); 
Noiseg = 0.01;
ng = sqrt(Noiseg);


zg = xg + (randn(1, Nt) + 1j .* rand(1, Nt))./sqrt(2) .* ng;

%% Signal avaialble
N_avail = length(t_avail);

x = exp(1j * (r0 + u * t_avail) * 4 * pi/lambda);

% Noise = sum(abs(x).^2)./(N_avail .* SNR); 
Noise = 0.01;
n = sqrt(Noise);

z = x + (randn(1, N_avail) + 1j .* rand(1, N_avail))./sqrt(2) .* n;

%% Doppler information

[~, ~, u_obs, sigma_obs, vel_axis, dv] = Spec(z(1:M), M, dt, lambda, SNR_db, 1, 1, 7);

%% Likelihood test

alpha0 = rand() * 5;
beta0 = rand() * 5;

% sigma_2u = betarnd(alpha0, beta0, 1);
sigma_2u = 1;
mu_0 = 10;


urange = linspace(-2,7, 5000);

for l = 1:length(urange)
    llu(l) = like(urange(l), z, t_avail, ng, mu_0, sigma_2u);
end

figure; plot(urange, llu);


%% MCMC algorithm

step_size = sigma_2u;
u_start = mu_0;


ll = like(u_start, z, t_avail, ng, mu_0, sigma_2u);
Nstep = 100000;
accept = [];
reject = [];
itern = [];

trailing_guardcells = 100;
Eps = 1;
epsilon = 2e-4;
sig(1) = sigma_2u;

w=[]; c = [];

for s = 1:Nstep
%     disp(s);
    u_new = get_proposal(u_start, step_size);
    u_new_debug(s) = u_new;
    ll_new = like(u_new, z, t_avail, ng, mu_0, sigma_2u);
    accept_r = exp(ll_new - ll);
    ran = rand;
    if ran < accept_r
        accept = [accept u_new];
        u_start = u_new;
        ll = ll_new;
        itern = [itern s]; 
        
%         sigma_2u = betarnd(alpha, beta, 1);
%     
%         mu_0 = u_start;
%     
%         step_size = sigma_2u;
        
    else
        reject = [reject u_new];
    end
    
    alpha = rand() * 5;
    beta = rand() * 5;

    
    %% Correlation check
    
    
    if s > trailing_guardcells
    
        [c, w] = change_sigma(u_new_debug, trailing_guardcells, epsilon, u_obs, sigma_obs, w, c);
        
        if length(c) >= 1
            if c(end) == 0

                disp(s);

                Eps = Eps+0.1;
            else
                break;
            end
        end

    end
    
    
    sigma_2u =  Eps * betarnd(alpha, beta, 1);
    
    mu_0 = u_start;
    
    step_size = sigma_2u;
    
    
    
end

figure; plot(reject); hold on; plot(accept);


burnin = round(0.25 * length(accept));
mu_ac = mean(accept(burnin:end)); 
%mu_re = mean(reject);

%% U debug

figure; plot(u_new_debug);

%% Correlation
var = abs(u_new_debug - mean(u_new_debug))./(norm(u_new_debug));
corrs = xcorr(var, var); 
corrs = corrs(2:end-1);
corrs = corrs(floor(length(corrs)/2):end);

figure; plot(corrs);

%% Plot the weights 

figure; plot(log(w));