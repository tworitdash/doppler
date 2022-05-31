%% HD signal generator

close all;
clear;

SNR_db = 20;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03; %0.03;              % Center frequency of the radar
c = 3e8;

x0 = 100;  % Start position of the target
a = 1;

u1 = 4;  u2 = 5;                    % Ground truth velocity of the target
a1 = 1; a2 = 4;

gt = [u1 u2 a1 a2];

NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
n_rot = 10;                % Number of rotations of radar

% Sections = [1 20 50 100 200];
% Sections = round(linspace(1, 200, 200));
Sections = 1;

% Sections = 1;
es = 0.5;

NMc = 2;



for l = 1:length(Sections)
for m = 1:NMc

disp('=========================================================================================');
    
disp(m);


disp('=========================================================================================');

NSec = Sections(l);                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
v_amb = lambda/(2*dT);
x1(1) = x0; % Initializing position
x2(1) = x0;

Z(1) = ( a1 .* exp(1j * 4 * pi/lambda .* x1(1)) + a2 .* exp(1j * 4 * pi/lambda .* x2(1)));

for i = 2:Nt
    x1(i) = x1(i - 1) + u1 * dT;
    x2(i) = x2(i - 1) + u2 * dT;
    Z(i) = ( a1 .* exp(1j * 4 * pi/lambda .* x1(i)) + a2 .* exp(1j * 4 * pi/lambda .* x2(i)) ) ;  % This for loop generates echo samples in time Nt times (HD)
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


% figure; plot(abs(Z_model)); hold on; plot(t_avail, real(Z_m), '*');

% figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

% This functions calculates the Doppler spectrum of the available samples
% in one rotation only N_Sweeps coherent samples. mu_obs and sigma_obs are
% the mean Doppler velocity and Doppler spectrum width with only N_Sweeps
% samples. Mu_obs is taken as the starting point for the MCMC algorithm.

% [~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);


%% MCMC parameters in a structure var

% E is the structure having options for MCMC

E.n = 4;                            % Number o fvariables to be estimated

% E.E0 = 2;%[mu_obs];                    % Initial value of u -> mu_obs
% E.sig = [3];                    % Initial value of the Std of the prior of u 

% E.H = [mu_obs + (sigma_obs-1.5)];
% E.L = [mu_obs - (sigma_obs - 1.5)];

E.L = [0 0 0 0]; 
E.H = [v_amb v_amb 5 5];

E.E0 = [2 2 1 1];

%% This section has MCMC algorithm 
% Num_Burnin = 5000;
% Num_MC     = 39000
Num_Simu = 10000000;
epsilon = 1e-2;
% Num_Simu = Num_Burnin + Num_MC; 

[accepted, rejected, itern, E_new, sample_MC_seq(l, m, :, :), sample_all(l, m, :, :), Flag] = MHu_uniformpriordouble(E, Num_Simu, Z_avail_vec, t_avail, x0, sigma_n, epsilon);


% Sample_mc_seq_with_burnin(l, m).seq = sample_MC_seq(l, :);

% 
% display(['Indices of accepted samples:'])
% itern
% 
% display(['Accepted samples:'])
% accepted
% 
% if Flag
%     sample_MC(l).seq = sample_MC_seq(l, NB(l)+1:end);
%     u_mean_unwraped(l, m) = mean(sample_MC(l).seq);
% %     u_mean(l, m) = mod(u_mean_unwraped(l, m), v_amb);
%     u_mean(l, m) = u_mean_unwraped(l, m);
%     if abs(u - u_mean(l, m)) < es
%         F(m) = 1;
%     else
%         F(m) = 0;
%     end
% end
% 
% end
% 
% FR(l) = sum(F)/NMc;
% 
% u_axis = linspace(-20, 20, 10000);
% 
% for o = 1:10000
%     ll(o) = LLu(u_axis(o), Z_avail_vec, t_avail, x0, sigma_n, 1);
% end

% figure(101); hold on; plot(u_axis, ll);
% 
end
sample_MC_seq_avg = squeeze(mean(sample_MC_seq, 2));
sample_all_avg = squeeze(mean(sample_all, 2));

    for k = 1:E.n
        figure(l+k); plot(sample_all_avg(k, :), 'Color', 'y');
        hold on; plot(sample_MC_seq_avg(k, :), 'Color', 'k');
        hold on; plot(ones(1, Num_Simu) .* gt(k),'--', 'Color', 'g');
    end

end



%%
% figure; 
% for l = 1:NMc
%     hold on; plot(Sample_mc_seq_with_burnin(1, l).seq);
% end
Nu = 200;
Na = 1000;

u_axis = linspace(0, v_amb, Nu);
a_axis = linspace(0.5, 1.5, Na);

for i = 1:Nu
    for k = 1:Na
        LL(i, k) = LLudouble([u_axis(i) a_axis(k)], Z_avail_vec, t_avail, x0, sigma_n, 1);
    end
end


figure; surface(u_axis, a_axis, LL.'); shading flat; colormap('jet');


