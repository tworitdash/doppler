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
n_rot = 10;                % Number of rotations of radar

% Sections = [1 20 50 100 200];
% Sections = round(linspace(1, 200, 200));


Sections = randi([20-2, 20+2], 1, n_rot);



% Sections = 1;
es = 0.5;

NMc = 2;

disp('=========================================================================================');
    
% disp(m);


disp('=========================================================================================');

NSec = [0 Sections];                  % Number of sections

Nt = NSweep * sum(NSec); % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
v_amb = lambda/(2*dT);

Gap = NSweep*(Sections - 1);
u_amb_gap = lambda./(2*Gap*dT);

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


 [~, ~, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dT, lambda, SNR_db, 1, 1, 7);


%% Available samples [Measurement model with only a few samples]

for k = 1:n_rot
    Z_model_re(k).data = Z_model(NSweep*sum(NSec(1:k)) + [1:NSweep*NSec(k+1)]);
    Z_avail_vec((k-1)*NSweep + [1:NSweep]) = Z_model_re(k).data(1:NSweep);
end

Z_avail = Z_model_re(1:NSweep, :);              % Takes only N_Sweep number of samples from all rotations             
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]); % vectorize all samples available

Z_avail_vec = Z_avail_vec_; %  + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [0:NSweep-1]; % This for loop calculates the time instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances

% [~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec(1:NSweep).', NSweep, dT, lambda, SNR_db, 1, 1, 5);