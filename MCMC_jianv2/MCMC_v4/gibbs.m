clear;
close all;

%% Gibbs sampler for 2 variables

v = 4;

lambda = 0.03;
f = v*2/lambda;

dt = 1e-3;

% T = 1/f;
nr= 100;
np = 10;
ns = 1;
t = eps:dt:(ns*np*nr)*dt;

% t = linspace(eps, 20*T, ns*np*nr);

% dt = t(2) - t(1);

fmax = 1/(2*dt);
for i = 1:nr
    t_avail(:, i) = t((i-1)*ns*np+1:(i-1)*ns*np+np);
end

t_a = reshape(t_avail, [1 nr*np]);

%% 
nt = length(t);
na = length(t_a);
SNR_db = 20;
SNR = 10^(SNR_db/10);
% x0 = 100;

s = cos(2 * pi * f * t);
sig_n2 = sum(abs(s).^2)/(nt * SNR);

z = s + sqrt(sig_n2) * (randn(1, nt)); % + 1j * rand(1, nt));

s_a = cos(2 * pi * f * t_a);
z_a = s_a + sqrt(sig_n2) * (randn(1, na)); % + 1j * rand(1, na));

% figure; plot( linspace(-fmax/2, fmax/2, nt), abs(fftshift(fft(z))).^2 );
% figure; plot( linspace(-fmax/2, fmax/2, na), abs(fftshift(fft(z_a))).^2 );
%%

% [~, ~, mu, sigma, vel_axis, dv] =  Spec2(z, nt, dt, lambda, SNR_db, 1, 1, 5);
% [~, ~, mu_a, sigma_a, vel_axis_a, dv_a] =  Spec2(z_a, na, dt, lambda, SNR_db, 1, 1, 8);

%% MCMC with Gibbs sampling: 
%% Task: Estimate f and sig_n from the data with MCMC Gibbs sampling 

K = 10000;
Ms = 10;
nu = na - 1;
f0 = 233.33;
omega0 = 2*pi*f0;
B10 = 1.1;
sig02 = sig_n2;


for m = 1:Ms
    disp(m);
    for k = 1:K
        D = z_a;
        B1khat = sum(D .* cos(omega0.*t_a))/sum((cos(omega0.*t_a)).^2);
        X1 = cos(omega0 .* t_a).';
        cov_B1 = sig02 .* inv(X1.' * X1);
        
        B1(k, m) = normrnd(B1khat, sqrt(cov_B1), 1);
        
        omega0hat = omega0;
        
        for l = 2:50
            J = -(B1(k, m) .* t_a .* sin(omega0hat(l-1) .* t_a)).';
            func = (B1(k, m) .* cos(omega0hat(l-1) .* t_a)).';
            delD = D.' - func;
            delomega = inv(J.' * J) * J.' * delD;
            omega0hat(l) = omega0hat(l-1)  + delomega;
        end
        
%         figure(100); hold on; plot(omega0hat/(2*pi));
        
        
        X_omega = -t_a.' .* B1(k, m) .* sin(omega0hat(end).*t_a).';
%         Dhat = B1(k) .* cos(omega0.*t_a);
%         s_omega2 = 1/(na - 1) .* (D - Dhat).' * (D - Dhat);
        cov_Omega = sig02 * inv(X_omega.' * X_omega);
        
        Omega(k, m) = normrnd(omega0hat(end), sqrt(cov_Omega), 1);
        
        omega0 = Omega(k, m);
        B10 = B1(k, m);
        freq(k, m) = mod(Omega(k, m)/(2*pi), fmax);
        vel(k, m) = mod(freq(k, m)*lambda/2, fmax*lambda/2);
    end
end


figure; plot(B1(:, 1)); hold on; plot(ones(1, K));
figure; plot(vel(:, 1)); hold on; plot(ones(1, K).*v)
figure; 
histogram(B1(:, 1)); 
figure; histogram(vel(:, 1));