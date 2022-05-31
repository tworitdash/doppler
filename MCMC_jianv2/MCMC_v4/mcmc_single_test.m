%% HD signal generator

close all;
clear;

SNR_db = 0;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03; %0.03;              % Center frequency of the radar
c = 3e8;

x0 = 100;                   % Start position of the target
u = 4;                      % Ground truth velocity of the target

NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
n_rot = 100;                % Number of rotations of radar

Sections = [1 20 50 100 200];
% Sections = round(linspace(1, 200, 200));

% Sections = 1;
es = 0.5;

NMc = 16;



for l = 1:length(Sections)
for m = 1:NMc

disp('=========================================================================================');
    
disp(m);


disp('=========================================================================================');

NSec = Sections(l);                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
v_amb = lambda/(2*dT);
x(1) = x0;                  % Initializing position 

Z(1) = exp(1j * 4 * pi/lambda .* x(1)); % First radar echo in time

for i = 2:Nt
    x(i) = x(i - 1) + u * dT;
    Z(i) = exp(1j * 4 * pi/lambda .* x(i));  % This for loop generates echo samples in time Nt times (HD)
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise/2);

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt)); % Adding complex noise 

% This function calculates the Doppler moments and plots the spectrum.
% Here, the ground truth signal is given as input
% mu is mean Doppler velocity, sigma is Doppler spectrum width, vel_axis is
% the velocity axis with Nt points. The last 3 inputs are for plot options

%  [~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available samples [Measurement model with only a few samples]


% Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]);
% Z_avail = Z_model_re(1:NSweep, :);
% Z_avail_vec = reshape(Z_avail, [NSweep * n_rot 1]);

% figure; histogram(real(Z_model)); figure; histogram(imag(Z_model));

Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]); 
Z_avail = Z_model_re(1:NSweep, :);              % Takes only N_Sweep number of samples from all rotations             
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]); % vectorize all samples available

% Noise_avail = sum(abs(Z_avail_vec_).^2)./(length(Z_avail_vec_) .* SNR);
% sigma_n_avail = sqrt(Noise_avail);

Z_avail_vec = Z_avail_vec_; %  + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [0:NSweep-1]; % This for loop calculates the time instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances


% figure; plot(t_avail, real(Z_avail_vec_)); hold on; plot(t_avail, real(Z_m), '*');

% figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

% This functions calculates the Doppler spectrum of the available samples
% in one rotation only N_Sweeps coherent samples. mu_obs and sigma_obs are
% the mean Doppler velocity and Doppler spectrum width with only N_Sweeps
% samples. Mu_obs is taken as the starting point for the MCMC algorithm.

% [~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);


%% MCMC parameters in a structure var

% E is the structure having options for MCMC

E.n = 1;                            % Number o fvariables to be estimated

E.E0 = 2;%[mu_obs];                    % Initial value of u -> mu_obs
E.sig = [3];                    % Initial value of the Std of the prior of u 

% E.H = [mu_obs + (sigma_obs-1.5)];
% E.L = [mu_obs - (sigma_obs - 1.5)];

%% This section has MCMC algorithm 
% Num_Burnin = 5000;
% Num_MC     = 39000
Num_Simu = 20000;
epsilon = 1e-2;
% Num_Simu = Num_Burnin + Num_MC; 
[accepted, rejected, itern, E_new, sample_MC_seq(l, :), sample_all, Flag, NB(l)] = MHu(E, Num_Simu, Z_avail_vec, t_avail, x0, sigma_n, epsilon);
% 
% display(['Indices of accepted samples:'])
% itern
% 
% display(['Accepted samples:'])
% accepted

if Flag
    sample_MC(l).seq = sample_MC_seq(l, NB(l)+1:end);
    u_mean_unwraped = mean(sample_MC(l).seq);
    u_mean(l, m) = mod(u_mean_unwraped, v_amb);
    if abs(u - u_mean(l, m)) < es
        F(m) = 1;
    else
        F(m) = 0;
    end
end

end

FR(l) = sum(F)/NMc;

u_axis = linspace(-20, 20, 10000);

for o = 1:10000
    ll(o) = LLu(u_axis(o), Z_avail_vec, t_avail, x0, sigma_n, 1);
end

figure(101); hold on; plot(u_axis, ll);

end

%%

