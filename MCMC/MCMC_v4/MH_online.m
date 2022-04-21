%% HD signal generator

close all;
clear;

SNR_db = 30;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03;             

r0 = 0;                     
u = 4;                      % Ground truth u

Nt = 10000;                 % Truth samples

          

dT = 1e-3;                  % t step
r(1) = r0;                  

Z(1) = exp(1j * 4 * pi/lambda .* r(1)); 

for i = 2:Nt
    r(i) = r(i - 1) + u * dT;
    Z(i) = exp(1j * 4 * pi/lambda .* r(i)); % Ground truth samples
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise);

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2); % Adding complex noise 


%% Available samples [Measurement model with only a few samples]

NSweep = 5;  


Z_model_re = reshape(Z, [100 100]); 
Z_avail = Z_model_re(1:NSweep, :); 

Z_avail_vec_ = reshape(Z_avail, [500 1]); % available samples for measurements


Z_avail_vec = Z_avail_vec_ + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:100
    t(:, k) = (k - 1) * 20 * NSweep + [1:NSweep]; % This for loop calculates the t instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances


%% MCMC parameters in a structure var

% E is the structure having options for MCMC

E.n = 1;                            % Number of variables to be estimated

E.E0 = [3.5354];                    % Initial value of u 
E.sig = [50000];                    % Initial value of the Std of the prior of u 


%% This section has MCMC algorithm 
No_iter = 10000; % Number of iterations

[accepted, rejected, itern, E_new] = MHu(E, No_iter, Z_avail_vec, t_avail, r0, sigma_n);

%% Plot MCMC outputs 


for i = 1:E.n

    
   figure(1000+i);plot(rejected(:, i)); hold on; plot(accepted(:, i)); % Accepted and rejected values
   
   burnin = round(0.25 * length(accepted(:, i)));                      % 25% of the data is taken as burnin
   figure(2000+i); histogram(accepted(burnin+1:end, i), 100);          % histogram of accepted
   burninrej = round(0.25 * length(rejected(:, i)));                   % Burnin for rejected
   figure(3000+i); histogram(rejected(burninrej+1:end, i), 100);       % Burnin for accepted
   mu_re = mean(accepted(burnin+1:end, i));                            % Mean of accepted 
   rej_re = mean(rejected(burnin+1:end, i));                           % Mean of rejected
end


