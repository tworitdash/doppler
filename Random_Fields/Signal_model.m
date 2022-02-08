function [Z_model, ZFFT, v_amb, vel_axis, dv, Vtmean, Vtspread, vr, PT, mu, sigma] = Signal_model(Nt, dt, lambda, x0, y0, z0, r0, D_min, D_max, N0, D0, N, u_mean, v_mean, u_sigma, v_sigma, th, ph, SNR_db)



markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;
v_amb = lambda/(4 * dt);


%% Simulate radar echoes


x(1, :) = x0;
y(1, :) = y0;
z(1, :) = z0;

% D = linspace(eps, 8e-3, 1000);

% dD = 1e-5;
% D = 0.6:dD:5.8; % in mm 
D = linspace(D_min, D_max, N);

dD = D(2) - D(1);
ND = N0 * exp(-3.67 * D./D0);
D_int = ND.*D.^6.*dD;


%%
% 
% u = 1 * ones(1, N);
% v = 1 * ones(1, N);
u = normrnd(u_mean, u_sigma, [1 N]);
v = normrnd(v_mean, v_sigma, [1 N]);
w = 0 * normrnd(0, 0, [1 N]);

Wt = 9.65 - 10.3 .* exp(-600 .* D .* 1e-3) - w;


%% Theoretical mean of Terminal Fall velocity

Nvt = Wt .* D_int .* ND .* dD;
Dvt = D_int .* ND .* dD;

Vtmean = sum(Nvt)./sum(Dvt);

Vtspread = sqrt(sum((Wt - Vtmean).^2 .* D_int .* ND .* dD)./sum(Dvt));

%%
txt = ['Terminal Fall Velocity'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
dtext = ['D_0 = ', num2str(D0), ' [mm]'];
xl = 'Diameter [mm]';

f = figure(109); hold on; f.Position = [10 10 1000 1000];
color = colors(2).c;
yl =  ['velocity [m/s]'];

marker = markers(1);
plott2(D, Wt, xl, yl, txt, 3, dtext, color, marker);
%%
vr = u .* cos(th) .* cos(ph) + v .* cos(th) .* sin(ph) + Wt .* sin(th);

Ref = sum(D_int);

A0 = (Ref).^(1/2);
Z(1) = sum(A0 .* exp(1j .* 4 * pi / lambda * r0));

%%
txt = ['Drop Size Distribution'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
dtext = ['Z_r = ', num2str(db(Ref)/2), ' [dB]']; 
xl = 'Diameter [mm]';

f = figure(110); hold on; f.Position = [10 10 1000 1000];
color = colors(2).c;
yl =  ['N(D) D^6 10^{-5} mm^6 m^{-3} mm^{-1}'];

marker = markers(1);
plott2(D, D_int, xl, yl, txt, 3, dtext, color, marker);
SNR = 10^(SNR_db/10);

Z(1) = sum(A0 .* exp(1j .* 4 * pi / lambda * r0));

for i = 2:Nt
    x(i, :) = x(i - 1, :) + u .* dt;
    y(i, :) = y(i - 1, :) + v .* dt;
    z(i, :) = z(i - 1, :) + Wt .* dt;
    
%     figure(101); hold on; plot(x(i, :), y(i, :), '+');
    
    r_(i, :) = sqrt(x(i, :).^2 + y(i, :).^2 + z(i, :).^2);
    A(i, :) = (Ref).^(1/2);
    Z(i) = sum(A(i, :) .* exp(1j .* 4 * pi / lambda .* r_(i, :)));
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);
sigma_n = sqrt(Noise);

Z_model = Z + sigma_n .* (randn(1, Nt));

[ZFFT, PT, mu, sigma, vel_axis, dv] = Spec(Z_model, Nt, dt, lambda, SNR_db, 1, 1, 5);
% 
% Z_re = PT./(dv .* Nt);
% ref_re =  db(Z_re/N)/2;
end
