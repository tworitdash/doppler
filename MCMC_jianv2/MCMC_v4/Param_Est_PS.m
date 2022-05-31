%% HD signal generator

close all;
clear;

SNR_db = 20;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03; %0.03;              % Center frequency of the radar
c = 3e8;

x0 = 100;  % Start position of the target
a = 1;

u1 = 3;  u2 = 3.2;                % Ground truth velocity of the target

gt = [u1 u2];

NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
n_rot = 100;                % Number of rotations of radar

% Sections = [1 20 50 100 200];
% Sections = round(linspace(1, 200, 200));
Sections = 1;



% Sections = 1;
es = 0.5;

NMc = 2;

disp('=========================================================================================');
    
% disp(m);


disp('=========================================================================================');

NSec = Sections;                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
v_amb = lambda/(2*dT);

Gap = NSweep*(Sections - 1) + 1;
u_amb_gap = lambda/(2*Gap*dT);

x1(1) = x0; % Initializing position
x2(1) = x0;

Z(1) = a * ( exp(1j * 4 * pi/lambda .* x1(1)) + exp(1j * 4 * pi/lambda .* x2(1)));

for i = 2:Nt
    x1(i) = x1(i - 1) + u1 * dT;
    x2(i) = x2(i - 1) + u2 * dT;
    Z(i) = a * ( exp(1j * 4 * pi/lambda .* x1(i)) + exp(1j * 4 * pi/lambda .* x2(i)) ) ;  % This for loop generates echo samples in time Nt times (HD)
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

%% Genetic Algorithm to solve the frequencies in a signal


Opt.fitnessfcn = @(x) Max_Likelihood(x(1:2), Z_avail_vec, t_avail, x0, sigma_n, 1);
% Opt.objective = @(x) Max_Likelihood(x(1:2), Z_avail_vec, t_avail, x0, sigma_n, 1);

% Opt.solver = 'patternsearch';
Opt.solver = 'ga';

lb = [0 0];
ub = [15, 15];
R0 = [2.7 3.5];

Opt.lb = [lb.'];
Opt.ub = [ub.'];
IP = [R0];
Opt.x0 = IP;
Opt.nvars = 2;
Opt.fitnesslimit = 50;

% Opt.options = optimoptions(@patternsearch,'PlotFcn', ...
%     {'psplotbestf', 'psplotbestx'},...
%      'UseParallel', true, 'MaxFunEval', Inf, 'MaxIter', 10000, 'MeshTolerance', 0, 'StepTolerance', 0);
% 
tic;
Opt.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', [IP], 'UseParallel',...
    true, 'FitnessLimit', Opt.fitnesslimit, 'FunctionTolerance', 0, 'MaxGenerations', 1000000);
[U, fval,exitFlag,output] = ga(Opt);

time_consumed = toc;
% tic;
% 
% [Result, fval, exf2, out] = patternsearch(Opt);
% time_consumed = toc;


% Opt.options = optimoptions('simulannealbnd','PlotFcns',...
%           {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
% 
% 
% [U, fval,exitFlag,output] = simulannealbnd(Opt.objective, Opt.x0, Opt.lb, Opt.ub, Opt.options);
% time_consumed = toc;



%% Likelihood test
Nu = 1000;
Nv = 1000;

u_test = linspace(0, 15, Nu);
v_test = linspace(0, 15, Nv);

parfor ui = 1:Nu
   for vi = 1:Nv

    LL(ui, vi) = Max_Likelihood([u_test(ui) v_test(vi)], Z_avail_vec, t_avail, x0, sigma_n, 1);
    
   end
end


txt = ['Number of Gap samples: ', num2str(Gap-1), ', Gap Ambiguity: ', num2str(u_amb_gap), ' [m/s]'];
figure; surface(u_test, v_test, -LL) %'LineWidth', 2, 'DisplayName', txt); 
shading flat; colormap('jet');
title(txt);

[minval, idx] = min(LL(:));
[row, column] = ind2sub(size(LL), idx);


