clear; 
close all;

%% signal model
v = 4;
lambda = 0.03;
f = 4*v/lambda; 
n_rot = 100;
M = 5;

% nsec = round(linspace(1, 200, 200));
nsec = 1;
nmc = 128;

for s = 1:length(nsec)
for m = 1:nmc
   
sec = nsec(s);

n = n_rot * sec * M;

% t = linspace(0, 100/f, n);
% dt = t(2) - t(1);
dt = 1e-3;
% t = eps:dt:(n)*dt; 
t = linspace(eps, (n-1)*dt, n);

SNR_db = 0;
SNR = 10^(SNR_db/10);



x = sin(2*pi*f*t);

Noise = sum(abs(x).^2)/(n .* SNR);
sigma_n = sqrt(Noise);
z = x + (sigma_n) * randn(1, n);



% figure; plot(t, z);

n_avail = M * n_rot;

for k = 1:n_rot
    t_a(:, k) = (k - 1) * sec * M + [0:M-1]; % This for loop calculates the time instances of the available samples
end

t_avail_in = reshape(t_a, [n_avail 1]);
t_avail = t_avail_in * dt;

z_avail = z(t_avail_in+1).';

% hold on; plot(t_avail, z_avail);


%% initial parameters

mu_start = 2*2/lambda;
sigma_start = 3*2/lambda;
f_start = normrnd(mu_start, sigma_start);


%% mcmc
n_iter = 20000;
burnin = 10000;


for i = 1:n_iter
    f_new = normrnd(f_start, sigma_start, 1);
    pf = p(f_start, f_start, sigma_start);
    pfnew = p(f_new, f_start, sigma_start);
    
    fl = likeli(f_start, z_avail, sigma_n, t_avail, pf);
    flnew = likeli(f_new, z_avail, sigma_n, t_avail, pfnew);
    
    if acceptance(fl, flnew)
        mc_seq(i, m) = f_new;
        all_seq(i, m) = f_new;
        f_start = f_new;
    else
        mc_seq(i, m) = f_start;
        all_seq(i, m) = f_new;
    end
end

figure(101); hold on; plot(1:n_iter, mc_seq(:, m)); grid on;
 
meanf(s, m) = mean(mc_seq(burnin+1:end, m));
if abs(meanf(s, m) - f) < 10
    flag(s, m) = 1; 
else
    flag(s, m) = 0;
end
end
    fratio(s) = sum(flag(s, :))/nmc;
end



% figure;plot(1:n_iter, all_seq, '*');hold on; plot(1:n_iter, mc_seq, '+');  hold on; plot(1:n_iter, f * ones(1, n_iter)); 

%% plot likelihood

f_axis = linspace(-1600, 1600, 200000);

for k = 1:200000
    llf(k) = likeli(f_axis(k), z_avail, sigma_n, t_avail, 1);
end

figure; plot(f_axis, llf);

%% likelihood

function ll = likeli(f, z, sigma_n, t_avail, pu)
    x = sin(2*pi*f*t_avail);
    nt = length(t_avail);
%     noise = z - x;
%     n_std = noise/sqrt(sigma_n);
    ll = log(pu) -nt * log(2*pi) - nt*2*log(sigma_n) - 1/(2*sigma_n^2) * (z - x).'  * (z - x) ;
%     llzu = -1/2 * log(2*pi) - 1/2 * n_std.^2;
%     ll = sum(llzu) + log(pu);
end

function pu = p(u, u_mean, u_sigma)
    pu = (normpdf(u, u_mean, u_sigma));
end

