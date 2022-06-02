%% HD signal generator

close all;
clear;

SNR_db = 0;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03; %0.03;              % Center frequency of the radar
c = 3e8;

x0 = 0;                   % Start position of the target
u = 4;                      % Ground truth velocity of the target

NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
% N_rot = round(linspace(1, 100, 10));
N_rot = 100;
            
for ni = 1:length(N_rot)

n_rot = N_rot(ni);                % Number of rotations of radar

% Sections = [1 20 50 100 200];
% Sections = round(linspace(1, 200, 20));
Sections = 200;

es = 0.5;

NMc = 16;

% if Z
if ni > 1

    clear Z;
    clear Z_avail;
    clear Z_model;
    clear Z_model_re;
    clear Z_avail_vec;
    clear Z_avail_vec_;
    clear t_avail;
    clear t;
    
end
    
% end




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
%  [~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available samples [Measurement model with only a few samples]

Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]); 
Z_avail = Z_model_re(1:NSweep, :);              % Takes only N_Sweep number of samples from all rotations             
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]); % vectorize all samples available

Z_avail_vec = Z_avail_vec_; %  + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [0:NSweep-1]; % This for loop calculates the time instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances


% [~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);


%% MCMC parameters in a structure var

% E is the structure having options for MCMC

E.n = 1;                            % Number o fvariables to be estimated

E.L = 0; 
E.H = v_amb;

E.E0 = 2;

%% This section has MCMC algorithm 

Num_Simu = 20000;
epsilon = 1e-2;

% [accepted, rejected, itern, E_new, sample_MC_seq(l, :), sample_all, Flag, NB(ni, l, m)] = MHu_uniformprior(E, Num_Simu, Z_avail_vec, t_avail, x0, sigma_n, epsilon);
[accepted, rejected, itern, E_new, sample_MC_seq(l, :), sample_all, Flag] = MHu_uniformprior(E, Num_Simu, Z_avail_vec, t_avail, x0, sigma_n, epsilon);

Sample_mc_seq_with_burnin(ni, l, m).seq = sample_MC_seq(l, :);


if Flag
    sample_MC(ni, l).seq = sample_MC_seq(l, NB(ni, l, m)+1:end);
    sample_mc_all(ni, l).seq = sample_all(l, NB(ni, l, m)+1:end);
    u_mean_unwraped(ni, l, m) = mean(sample_MC(ni, l).seq);
    u_mean(ni, l, m) = u_mean_unwraped(ni, l, m); % mod(u_mean_unwraped(ni, l, m), v_amb);
    if abs(u - u_mean(ni, l, m)) < es
        F(ni, l, m) = 1;
    else
        F(ni, l, m) = 0;
    end
else
    NB(ni, l, m) = 2000;
    sample_MC(ni, l).seq = sample_MC_seq(l, NB(ni, l, m)+1:end);
    sample_mc_all(ni, l).seq = sample_all(l, NB(ni, l, m)+1:end);
    u_mean_unwraped(ni, l, m) = mean(sample_MC(ni, l).seq);
    u_mean(ni, l, m) = u_mean_unwraped(ni, l, m); %mod(u_mean_unwraped(ni, l, m), v_amb);
    F(ni, l, m) = 0;
end

if m == 1
    figure; plot(sample_mc_all(ni, l).seq); hold on; plot(sample_MC(ni, l).seq); hold on; ...
    plot(u .* ones(1, length(sample_MC(ni, l).seq)), '--', 'Color', 'g');
end

end

FR(ni, l) = sum(F(ni, l, :), 3)/NMc;

RMS(ni, l) = sqrt( sum((squeeze(u_mean(ni, l, :)) - u).^2)/NMc ); 

end

end

%% Plot the error


xl = 'Number of scans [N_{scans}]';
yl = 'Number of echo samples in gaps [N_{sec} - N_[pulse]]';
zl = 'Root mean square error [RMS]';
tl = ['RMS error with ', num2str(NMc), ' Monte Carlo Runs'];
surplot_pcolor(N_rot, NSweep*(Sections - 1), RMS.', xl, yl, zl, tl);

