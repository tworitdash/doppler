%% HD signal generator

close all;
clear;

SNR_db = 30;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 0.03;              % Center frequency of the radar

x0 = 0;                   % Start position of the target
u = 4;                      % Ground truth velocity of the target

NSweep = 5;                 % Number of sweeps available per beamwidth 
                            %(low resolution, used for measurement)
n_rot = 10;                % Number of rotations of radar
NSec = 200;                  % Number of sections

Nt = n_rot * NSweep * NSec; % Number of samples in the ground truth (HD signal)


dT = 1e-3;                  % PRT of radar
vamb = lambda/(2*dT);

x(1) = x0;                  % Initializing position 

Z(1) = exp(1j * 4 * pi/lambda .* x(1)); % First radar echo in time

for i = 2:Nt
    x(i) = x(i - 1) + u * dT;
    Z(i) = exp(1j * 4 * pi/lambda .* x(i));  % This for loop generates echo samples in time Nt times (HD)
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             % noise variance with the data and given SNR
sigma_n = sqrt(Noise);

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2); % Adding complex noise 

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

Z_model_re = reshape(Z, [NSweep * NSec n_rot]); 
Z_avail = Z_model_re(1:NSweep, :);              % Takes only N_Sweep number of samples from all rotations             
Z_avail_vec_ = reshape(Z_avail, [NSweep * n_rot 1]); % vectorize all samples available

% Noise_avail = sum(abs(Z_avail_vec_).^2)./(length(Z_avail_vec_) .* SNR);
% sigma_n_avail = sqrt(Noise_avail);

Z_avail_vec = Z_avail_vec_ + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:n_rot
    t(:, k) = (k - 1) * NSec * NSweep + [1:NSweep]; % This for loop calculates the time instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances


% figure; plot(t_avail, real(Z_avail_vec_)); hold on; plot(t_avail, real(Z_m), '*');

% figure; plot(t_avail, real(Z_avail_vec)); hold on; plot(t_avail, imag(Z_avail_vec))

% This functions calculates the Doppler spectrum of the available samples
% in one rotation only N_Sweeps coherent samples. mu_obs and sigma_obs are
% the mean Doppler velocity and Doppler spectrum width with only N_Sweeps
% samples. Mu_obs is taken as the starting point for the MCMC algorithm.

na = length(Z_avail_vec);

[~, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_avail_vec.', na, dT, lambda, SNR_db, 1, 1, 5);


%% Deterministic Annealing

E.n = 1; % Numnber of parameters to be optimized

E.L = 0;
E.H = vamb;


%% Create clusters randomly using principles of random fields
% nuaxis = 1000;
% 
% u_axis = linspace(E.L, E.H, nuaxis);
% dx = u_axis(2) - u_axis(1);
% 
% f = 0.02;
% lambdax = 1;
% kx0 = 2*pi/lambdax;
% nkx = 100;
% kx = linspace(-kx0, kx0, nkx);
% 
% dkx = kx(2) - kx(1);
% 
% nx1 = 1 + floor(0.9999.*f*nkx);
% nx2 = nkx - nx1 + 1;
% 
% F = (rand(1, nuaxis));
% 
% for i = 1:length(kx)
%     F_(i) = sum(F .* exp(-1j .* (kx(i) .* u_axis)) .* dkx);
% end
% 
% F_(nx1:nx2) = 0;
% 
% for i = 1:length(u_axis)
%     F__(i) = sum(F_ .* exp(1j .* (u_axis(i) .* kx)) .* dx);
% end
% 
% 
% Z_zeropad = Z_avail_vec;
% Z_zeropad(end+1:nuaxis) = 0;
% 
% 
% [ZFFT, ~, mu_obs, sigma_obs, vel_axis_obs, dv_obs] =  Spec(Z_zeropad.', nuaxis, dT, lambda, SNR_db, 1, 1, 8);
% 
% figure; plot(u_axis, db(abs(F__)./(max(abs(F__)).*ZFFT)), 'LineWidth', 2); grid on;

%% making clusters for the parameter of estimation

nC = 10;
nCN = 1000; 
vChunks = (E.H - E.L)/10;
vChunkAxis = E.L:vChunks:E.H;

for k = 2:length(vChunkAxis)
    vClusters(k-1, :) = linspace(vChunkAxis(k - 1), vChunkAxis(k), nCN);
end


beta = 1;
Niter = 100000;

for iter = 2:Niter

    for ti = 1:na
        for l = 1:nC
            for m = 1:nCN
                y(l, m, ti) = exp( 1j .* 4*pi/lambda .*( x0 + vClusters(l, m) .* t_avail(ti) ) );
                zri = [real(Z_avail_vec(ti)) imag(Z_avail_vec(ti))];
                yri = [real(y(l, m, ti)) imag(y(l, m, ti))];
                E_lmt(l, m, ti) = [zri - yri] * [zri - yri].';
                PrzC_num(l, m, ti) = exp(-beta.*abs(E_lmt(l, m, ti)));
            end
        end

        PrzC_den(ti) = sum(sum(squeeze(PrzC_num(:, :, ti))));
        PrzC(:, :, ti) = squeeze(PrzC_num(:, :, ti))./PrzC_den(ti);
    end

end




