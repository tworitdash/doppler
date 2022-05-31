clear; 
close all;

lambda = 0.03; 
dt = 1e-3;

u = 3;
SNR_db = 30;
SNR = 10^(SNR_db/10);

n_rot = 10;
K = linspace(0, n_rot-1, n_rot);

M = 5;
Nsec = 2;
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

%% Likelihood test
sigma_2u = 1;
mu_0 = 0;


urange = linspace(-2,7, 5000);

for l = 1:length(urange)
    llu(l) = like(urange(l), z, t_avail, ng, mu_0, sigma_2u);
end

figure; plot(urange, llu);


%% MCMC algorithm



step_size = 0.0001;

u_start = 2;


ll = like(u_start, z, t_avail, ng, mu_0, sigma_2u);
Nstep = 100000;
accept = [];
reject = [];

for s = 1:Nstep
    u_new = get_proposal(u_start, step_size);
    ll_new = like(u_new, z, t_avail, ng, mu_0, sigma_2u);
    accept_r = exp(ll_new - ll);
    ran = rand;
    if ran < accept_r
        accept = [accept u_new];
        u_start = u_new;
        ll = ll_new;
    else
        reject = [reject u_new];
    end
end

figure; plot(reject); hold on; plot(accept);


burnin = round(0.25 * length(accept));
mu_ac = mean(accept(burnin:end)); 
%mu_re = mean(reject);



