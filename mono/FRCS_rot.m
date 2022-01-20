%% Fast rotation convergence study

clear;
% close all;

m_ = load('markers.mat'); 
markers = m_.markers;


dr = 2.27;
r = eps:dr:1e3;
Nr = length(r);
R = 2.27;


r_avg_op = 0;
multi_rot_op = 1;


Omega_rpm = 60;
n_rot = 200;
Nrot = 1:n_rot;
PRT = 1e-3;
Vxm = 3 .* ones(size(r));
Vym = 0 .* ones(size(r));

Vxs = 1.5 .* ones(size(r));
Vys = 0 .* ones(size(r));

phi_0_deg = 0;
theta = 0;
n_sig = sqrt(0.001);
lambda = 0.03;
BW_deg = 1.8;

parfor i = 1:n_rot

[s, v_amb, dv, SNR, hs, N, t, mu_prior, phi_axis, Nphi, sec] = TDG(r, Omega_rpm, Nrot(i), lambda, BW_deg, PRT, Vxm, Vxs, Vym, Vys, phi_0_deg, theta, n_sig);


[mu, sigma, mu_r_sigma, sigma_r_sigma, r_new, dv] = retrieval(s, v_amb, hs, Nr, Nphi, sec, r_avg_op, r, R, multi_rot_op, mu_prior);

mu_phiwind(i) = mu(1, 1);
sigma_phiwind(i) = sigma(1, 1);
SNR_avg(i) = SNR(1, 1); 

if r_avg_op == 1
    r_new = r_new(1:end-1);
end

end

%% True Mean velocity at the direction of wind 

mu_true = Vxm(1) .* ones(1, n_rot);
sigma_true = Vxs(1) .* ones(1, n_rot);

%% Plot for convergence with range averaging

txt = [' \Omega = ', num2str(Omega_rpm), ' [rpm]'];

dtext = ['\mu_{retrieved} with range averaging'];
xl = 'Number of range cells averaged';

f = figure(105); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['\mu_{retrieved} at \phi = \phi_{wind}'];
plott(Nrot, mu_phiwind, xl, yl, txt, 2, dtext, color, markers(1))

dtext = ['\mu_{True}'];

hold on; f.Position = [10 10 1000 1000];

plott(Nrot, mu_true, xl, yl, txt, 2, dtext, color, markers(2))

txt = [' \Omega = ', num2str(Omega_rpm), ' [rpm]'];

dtext = ['\sigma_{retrieved} with range averaging'];
xl = 'Number of range cells averaged';

f = figure(106); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['\sigma_{retrieved} at \phi = \phi_{wind}'];
plott(Nrot, sigma_phiwind, xl, yl, txt, 2, dtext, color, markers(1))

dtext = ['\sigma_{True}'];

hold on; f.Position = [10 10 1000 1000];

plott(Nrot, sigma_true, xl, yl, txt, 2, dtext, color, markers(2))