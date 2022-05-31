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

for l = 1:length(Sections)
    

NSec = Sections(l);                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
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

 [~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


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

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);


%% MCMC parameters in a structure var

% E is the structure having options for MCMC

E.n = 1;                            % Number o fvariables to be estimated

E.E0 = 2;%[mu_obs];                    % Initial value of u -> mu_obs
E.sig = [3];                    % Initial value of the Std of the prior of u 

% E.H = [mu_obs + (sigma_obs-1.5)];
% E.L = [mu_obs - (sigma_obs - 1.5)];

%% This section has MCMC algorithm 
Num_Burnin = 1000;
Num_MC     = 19000
Num_Simu = Num_Burnin + Num_MC; 
[accepted, rejected, itern, E_new, sample_MC_seq(l, :), sample_all] = MHu(E, Num_Simu, Z_avail_vec, t_avail, x0, sigma_n);

display(['Indices of accepted samples:'])
itern

display(['Accepted samples:'])
accepted

% close all
%% Plot MCMC outputs 


figure(101); hold on;

plot(1:Num_Simu,sample_MC_seq(l, :), 'linewidth', 2, 'DisplayName', ['N_{s} - N_{pulse} =  ', num2str(NSec*NSweep - NSweep)]);

    
% lgd = legend('Monte Carlo samples','MC sequence','Gound truth');


xlabel('Indices of all the samples',  'FontSize', 30)
ylabel('Value of the sequences',  'FontSize', 30)
title('Sequences of MC',  'FontSize', 30)
grid on; set(gca, 'FontSize', 30, 'LineWidth', 4);
box on;

end

%% Histogram plots

nbins = linspace(0, 20, 20);

[h1, edges] = histcounts(sample_MC_seq(1, :), nbins);
centers = mean([edges(1:end-1); edges(2:end)]);
h2 = histcounts(sample_MC_seq(2, :), nbins);
h3 = histcounts(sample_MC_seq(3, :), nbins);
h4 = histcounts(sample_MC_seq(4, :), nbins);
h5 = histcounts(sample_MC_seq(5, :), nbins);

h_ground_truth = histcounts(u, nbins);


h1 = h1 / max(h1);
h2 = h2 / max(h2);
h3 = h3 / max(h3);
h4 = h4 / max(h4);
h5 = h5 / max(h5);


h_ground_truth_norm = h_ground_truth / max(h_ground_truth);

figure;

bar_handle = bar(centers', [h1', h2', h3', h4', h5', h_ground_truth_norm']);
bar_handle(1).FaceColor = 'b';
bar_handle(2).FaceColor = 'r';
bar_handle(3).FaceColor = 'g';
bar_handle(4).FaceColor = 'y';
bar_handle(5).FaceColor = 'm';
bar_handle(6).FaceColor = 'k';

lgd = legend(['N_{s} - M =  ', num2str(Sections(1)*NSweep - NSweep)], ...
            ['N_{s} - M =  ', num2str(Sections(2)*NSweep - NSweep)], ...
            ['N_{s} - M =  ', num2str(Sections(3)*NSweep - NSweep)], ...
            ['N_{s} - M =  ', num2str(Sections(4)*NSweep - NSweep)], ...
            ['N_{s} - M =  ', num2str(Sections(5)*NSweep - NSweep)],...
            'Ground truth');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
xlabel('u',  'FontSize', 30)
ylabel('',  'FontSize', 30)
title('Histogram of MC sequence',  'FontSize', 30)
grid on; set(gca, 'FontSize', 30, 'LineWidth', 4);
box on;



figure(101); 

hold on; 
plot(1:Num_Simu, u*ones(1,Num_Simu), 'g-.','linewidth',2, 'DisplayName', 'Grounnd Truth');
lgd = legend;
lgd.FontSize = 20;
lgd.FontWeight = 'bold';