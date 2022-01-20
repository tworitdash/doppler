%% Fast rotation convergence study

clear;
% close all;

m_ = load('markers.mat'); 
markers = m_.markers;


dr = 2.27;
r = eps:dr:100;
Nr = length(r);
N_avg = 1;
R = N_avg .* dr;
r_avg_op = 0;
multi_rot_op = 1;


Omega_rpm = 60;
n_rot = 16;
n_bw = 1;
PRT = n_bw * 1e-3;
Vxm = 0 .* ones(size(r));
Vym = 0 .* ones(size(r));

Vspread = linspace(eps, 7.5, 100);
% Vspread = [1 2];

for k = 1:length(Vspread)

Vxs = Vspread(k) .* ones(size(r));
Vys = 0 .* ones(size(r));

phi_0_deg = 0;
theta = 0;
n_sig = sqrt(0.001);
lambda = 0.03;
BW_deg = n_bw.*1.8;

for i = 1:length(N_avg)

[s, v_amb, dv, SNR, hs, N, t, mu_prior, phi_axis, Nphi, sec] = TDG(r, Omega_rpm, n_rot, lambda, BW_deg, PRT, Vxm, Vxs, Vym, Vys, phi_0_deg, theta, n_sig);


[mu(k, :, :), sigma(k, :, :), mu_r_sigma, sigma_r_sigma, r_new, dv] = retrieval_V2(s, v_amb, hs, Nr, Nphi, sec, r_avg_op, r, R(i), multi_rot_op, mu_prior, phi_axis);

mu_phiwind(i, k) = mu(k, 1, 1);
sigma_phiwind(i, k) = sigma(k, 1, 1);

if r_avg_op == 1
    r_new = r_new(1:end-1);
end

end
end
%% Plot for convergence with range averaging

txt = [' \Omega = ', num2str(Omega_rpm), ' [rpm]'];

dtext = ['\mu_{retrieved} in a Beam = ', num2str(n_bw), ' \times BW'];
xl = 'True spectrum width [m.sec^{-1}]';

f = figure(105); hold on; f.Position = [10 10 1000 1000];
color = 'b';
yl =  ['\mu_{retrieved} at \phi = \phi_{wind}'];
plott(Vspread, mu_phiwind(1, :), xl, yl, txt, 2, dtext, color, markers(1))


txt = [' \Omega = ', num2str(Omega_rpm), ' [rpm]'];

dtext = ['\sigma_{retrieved} in a Beam = ', num2str(n_bw), ' \times BW'];
xl = 'True spectrum width [m.sec^{-1}]';

f = figure(106); hold on; f.Position = [10 10 1000 1000];
color = 'b';
yl =  ['\sigma_{retrieved} at \phi = \phi_{wind}'];

plott(Vspread, sigma_phiwind(1, :), xl, yl, txt, 2, dtext, color, markers(1))


