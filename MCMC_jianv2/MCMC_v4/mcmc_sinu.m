clear; 
close all;

%% signal model
v = 3;
lambda = 0.03;
f = 2*v/lambda; 
n_rot = 100;
m = 5;
sec = 1;

n = n_rot * sec * m;

% t = linspace(0, 100/f, n);
% dt = t(2) - t(1);
dt = 1e-3;
% t = eps:dt:(n)*dt; 
t = linspace(eps, (n-1)*dt, n);
sigma_n = 0.01;

x = sin(2*pi*f*t);
z = x + sqrt(sigma_n) * randn(1, n);

figure; plot(t, z);

n_avail = m * n_rot;

for k = 1:n_rot
    t_a(:, k) = (k - 1) * sec * m + [0:m-1]; % This for loop calculates the time instances of the available samples
end

t_avail_in = reshape(t_a, [n_avail 1]);
t_avail = t_avail_in * dt;

z_avail = z(t_avail_in+1).';

hold on; plot(t_avail, z_avail);


%% initial parameters

mu_start = 300;
sigma_start = 100;
f_start = normrnd(mu_start, sigma_start);


%% mcmc
n_iter = 200000;



for i = 1:n_iter
    f_new = normrnd(f_start, sigma_start, 1);
    pf = p(f_start, f_start, sigma_start);
    pfnew = p(f_new, f_start, sigma_start);
    
    fl = likeli(f_start, z_avail, sigma_n, t_avail, pf);
    flnew = likeli(f_new, z_avail, sigma_n, t_avail, pfnew);
    
    if acceptance(fl, flnew)
        mc_seq(i) = f_new;
        all_seq(i) = f_new;
        f_start = f_new;
    else
        mc_seq(i) = f_start;
        all_seq(i) = f_new;
    end
end

figure;plot(1:n_iter, all_seq, '*');hold on; plot(1:n_iter, mc_seq, '+');  hold on; plot(1:n_iter, f * ones(1, n_iter)); 

%% plot likelihood

f_axis = linspace(-400, 400, 200000);

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

