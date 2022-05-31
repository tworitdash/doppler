clear;
close all;

dt = 0.001; 
M = 5; 
N = 100;

u_true = 2;
lambda = 0.03;
Noise = 0.01;
sigma_n = sqrt(Noise);
t = [];

for k = 0:9
    for i = 0:M-1
        t = [t N*k*dt + i*dt];
    end
end

Nlen = length(t);

x = exp(1j .* u_true * t);
n = randn(1, Nlen) + 1j .* randn(1, Nlen);

z = x + n;

sigma_2u = 1;
mu_0 = 0;

step_size = 0.001;

u_start = 1;


ll = like(u_start,  z, t, sigma_n, mu_0, sigma_2u);
Nstep = 100000;
accept = [];
reject = [];

for s = 1:Nstep
    u_new = get_proposal(u_start, step_size);
    ll_new = like(u_new, z, t, sigma_n, mu_0, sigma_2u);
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


urange = linspace(-2, 7, 5000);

for l = 1:length(urange)
    llu(l) = like(urange(l), z, t, sigma_n, mu_0, sigma_2u);
end

figure; plot(urange, llu);


function u = get_proposal(u, step_size)
    u = u + randn()*step_size;
end



function ll = like(u, z, t, sigma_n, mu_0, sigma_2u)
    x = exp(1j .* u * t);
    n = z - x;
    n_re_std = real(n) * sqrt(2) /sigma_n;
    n_im_std = imag(n) * sqrt(2) /sigma_n;
    
    ll1 = -1/2*log(2*pi) - 1/2 * n_re_std.^2;
    ll2 = -1/2*log(2*pi) - 1/2 * n_im_std.^2;
    
    log_pu = -1/2*log(2*pi*sigma_2u)-1/2*(u-mu_0)^2/sigma_2u;
    
    ll = sum(ll1) + sum(ll2)  + log_pu;
end


