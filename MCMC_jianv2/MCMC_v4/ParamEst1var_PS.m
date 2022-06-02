%% HD signal generator

close all;
clear;

SNR_db = 20;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03; %0.03;              % Center frequency of the radar
c = 3e8;

x0 = 0;  % Start position of the target
a = 1;

Nscatters = 1; 

u1mean = 4; 
u1spread = 1;


gt = [u1mean u1spread];

var(1).name = 'Mean Velocity';
var(2).name = 'Spectrum Width';

NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
n_rot = 100;                % Number of rotations of radar

% Sections = [1 20 50 100 200];
% Sections = round(linspace(1, 200, 200));
Sections = 200;



% Sections = 1;
es = 0.5;

% NMc = 2;

disp('=========================================================================================');
    
% disp(m);


disp('=========================================================================================');

NSec = Sections;                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
v_amb = lambda/(2*dT);

Gap = NSweep*(Sections - 1) + 1;
u_amb_gap = lambda/(2*Gap*dT);

% u = normrnd(u1mean, u1spread, [1 Nscatters]);  %u2 = 3.2;                % Ground truth velocity of the target
% u = 2 + 3 .* rand(1, Nscatters);
% u = u1mean - u1spread/2 + u1spread .* rand(1, Nscatters);
u = u1mean;

% x1(1, :) = x0*rand(1, Nscatters); % Initializing position
% x2(1) = x0;
x1(1) = x0;

Z(1) = sum( exp(1j * 4 * pi/lambda .* x1(1)));

for i = 2:Nt
    x1(i) = x1(i - 1) + u * dT;
%     x2(i) = x2(i - 1) + u2 * dT;
    Z(i) = sum( exp(1j * 4 * pi/lambda .* x1(i))) ;  % This for loop generates echo samples in time Nt times (HD)
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise/2);

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt)); % Adding complex noise 


[~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available samples [Measurement model with only a few samples]


Z_model_re = reshape(Z_model, [NSweep * NSec n_rot]); 
Z_avail = Z_model_re(1:NSweep, :);              % Takes only N_Sweep number of samples from all rotations             
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]); % vectorize all samples available

Z_avail_vec = Z_avail_vec_; %  + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [0:NSweep-1]; % This for loop calculates the time instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec.', length(Z_avail_vec), dT, lambda, SNR_db, 1, 1, 5);

%% Genetic Algorithm to solve the frequencies in a signal

NMc = 16;

for k = 1:NMc

Nscatters_p = 1;

umean_start = 2;
uspread_start = 2;

% ustart = normrnd(umean_start, uspread_start, [1 Nscatters_p]);
% ustart = (mu_obs - sigma_obs) + 2 * sigma_obs .* rand(1, Nscatters_p);
ustart = [umean_start]; % uspread_start];

% Opt.fitnessfcn = @(x) Max_LikelihoodGauss(x(1:2), Nscatters_p, Z_avail_vec, t_avail, x0, sigma_n, 1);
Opt.objective = @(x) Max_LikelihoodGaussvar1(x, Nscatters_p, Z_avail_vec, t_avail, x0, sigma_n, 1);

Opt.solver = 'patternsearch';
% Opt.solver = 'ga';

% lb = (mu_obs - sigma_obs) .* ones(1, Nscatters_p); % 0 .* ones(1, Nscatters_p);
% ub = (mu_obs + sigma_obs) .* ones(1, Nscatters_p); %15 .* ones(1, Nscatters_p);

lb = [mu_obs - sigma_obs]; % 0.5];
ub = [mu_obs + sigma_obs]; % 2.5];

R0 = ustart;

% lb = [0 0];
% ub = [15, 15];
% R0 = [2.7 3.5];

Opt.lb = [lb.'];
Opt.ub = [ub.'];
IP = [R0];
Opt.x0 = IP;
Opt.nvars = 1;
Opt.fitnesslimit = 50;

Opt.options = optimoptions(@patternsearch,'PlotFcn', ...
    {'psplotbestf', 'psplotbestx'},...
     'UseParallel', true, 'MaxFunEval', Inf, 'MaxIter', 100, 'MeshTolerance', 0, 'StepTolerance', 0);

% tic;
% Opt.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
%     'InitialPopulationMatrix', [IP], 'UseParallel',...
%     true, 'FitnessLimit', Opt.fitnesslimit, 'FunctionTolerance', 0, 'MaxGenerations', 1000);
% [Result, fval,exitFlag,output] = ga(Opt);
% 
% time_consumed = toc;
tic;
% 
[Result(k, :), fval(k), exf2(k), out(k)] = patternsearch(Opt);
time_consumed(k) = toc;


end


%% Error Calculation 


Nvar = Opt.nvars;

for er = 1:Nvar
    Error(er) = sqrt( sum((Result(er, :) - gt(er)).^2)/ NMc );
end




% Opt.options = optimoptions('simulannealbnd','PlotFcns',...
%           {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
% Opt.options.FunctionTolerance = 0;
% 
% 
% [Result, fval,exitFlag,output] = simulannealbnd(Opt.objective, Opt.x0, Opt.lb, Opt.ub, Opt.options);
% time_consumed = toc;

%% Histograms

for er = 1:Nvar
    figure; 
%     nbins = linspace(min(Result(:, er)), max(Result(:, er)), NMc/4);
    nbins = linspace(mu_obs-sigma_obs, mu_obs+sigma_obs, NMc);
    [n1, edges] = histcounts(Result(:, er), nbins);
    centers = mean([edges(1:end-1); edges(2:end)]);
    n2 = histcounts(gt(er), nbins);
    
    n1 = n1/max(n1);
    n2 = n2/max(n2);
    
%     hold on; histogram(gt(er), 50);figure;
%     figure;
    bar_handle = bar(centers', [n1',n2']);
    set(bar_handle, {'DisplayName'}, {'PS algorithm'; 'Ground Truth'});
    legend;
    title(['Histogram of ', var(er).name], 'FontSize', 16);
    xlabel(var(er).name, 'FontSize', 16);
    grid on;
    
    print(['PS_Histogram_var1_', num2str(er)], '-depsc');
    
end



%% Reconstruction


% u = 2 + 3 .* rand(1, Nscatters);

x1r(1) = x0; %*rand(1, Nscatters_p); % Initializing position
% x2(1) = x0;
% ur = mean(Result(:, 1)) - mean(Result(:, 2))/2 + mean(Result(:, 2)) .* rand(1, Nscatters_p);
ur = mean(Result);

Zr(1) = sum( exp(1j * 4 * pi/lambda .* x1r(1)));

for i = 2:Nt
    x1r(i, :) = x1r(i - 1, :) + ur * dT;
%     x2(i) = x2(i - 1) + u2 * dT;
    Zr(i) = sum( exp(1j * 4 * pi/lambda .* x1r(i))) ;  % This for loop generates echo samples in time Nt times (HD)
end

Noise = sum(abs(Zr).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise/2);

Z_modelr = Zr + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt)); % Adding complex noise 


[~, ~, mur, sigmar, vel_axisr, dvr] = Spec(Z_modelr, Nt, dT, lambda, SNR_db, 1, 1, 8);
print('PS_Doppler_Spectrum_var1', '-depsc');