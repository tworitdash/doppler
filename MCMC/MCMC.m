clear; 
close all;

%% Random data generator from a Normal distribution:

x = normrnd(10, 1, [1 3000]);
% figure; histogram(x);

obs = x(randi(3000, 1, 15));

figure; histogram(obs);

mu_obs = mean(obs);
std_obs = std(obs(1:5));
%% The model for MCMC

[accepted, rejected] = metropolis_hastings([mu_obs, 0.1], 50000, obs);

figure; plot(rejected); hold on; plot(accepted);

burnin_num = round(length(accepted) .* 0.25);


h = mean(accepted(burnin_num+1:end));

%% Comparison

ret = normrnd(mu_obs, h, [1 3000]);

figure; histogram(x); 
hold on; histogram(obs); 
hold on; histogram(ret);

% figure; histogram(accepted); 