clear;
close all;

m_ = load('markers.mat'); 
markers = m_.markers;
%% Generating time domain data for a monochriomatic wind within a rotating radar
dr = 2000;
R = 15e3;
% r = 0.1:dr:R;

r = 2e2;
Nr = length(r);

BW_deg = 1.8;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = linspace(1, 60, 4);
% Omega_rpm = 1;
lambda = 0.03;
BW = BW_deg * pi/180;
Omega = Omega_rpm * 2*pi/60; % rotation speed in rad/s
for m = 1:length(Omega_rpm)
    
%% Monte Carlo loop
MC = 32;


    
clear s; 

sec = ceil(2*pi/BW);

ph = 0:BW:(sec)*BW;

T = ph(end)/Omega(m);


Nf = floor(T/PRT);

hs = floor(Nf/sec);

N = hs*sec;

t = linspace(0, (Nf-1)*PRT, N);

% t = 0:PRT:(N - 1)*PRT;


hs_(m) = hs;

phi_0 = phi_0_deg * pi/180;
phi_axis_1 = phi_0:BW:phi_0+Omega(m)*t(end);
% phi_axis = linspace(phi_0, phi_0+Omega(m)*t(end), sec);

phi_axis = mean([ph(1:end-1); ph(2:end)]);
Nphi = length(phi_axis);
% phi_axis = zeros(size(phi_axis));

%% variables for signal model

beta_wind_deg = 0;
theta = 0;

beta_wind = beta_wind_deg .* pi/180;

v_amb = lambda/(4 * PRT);

Vmean = linspace(-v_amb, v_amb, 100);
% Vmean = 5;
% Vmean = 0;
Vspread = linspace(eps, 3, 100);
u_phi = 0; v_phi = 0;

s = zeros(MC, Nr, N);

for n = 1:length(Vmean)
    
for c = 1:MC
    
    
for i = 1:Nr
    V = normrnd(Vmean(n), 1, [1 N]);
    u = V .* cos(beta_wind);
    v = V .* sin(beta_wind);
    
    mu = V .* cos((beta_wind - Omega(m) .* t - phi_0)) .* cos(theta); % +  v .* sin(beta_wind - Omega .* t - phi_0) .* cos(theta);
    
%     p = Omega(m) .* t *180/pi;
%     figure(103); hold on; plot(p, mu, 'DisplayName', [num2str(Omega_rpm(m))]); legend;
    n_sig = sqrt(0.001); 
  
    [s(c, i, :), SNR(i)] = TD_generator(mu, lambda, beta_wind, phi_0, Omega(m), t, n_sig);
end

%% Retrieval

if (mod(hs, 2) == 0)
    hsnew = hs+1;
%     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
else
    hsnew = hs;
    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
end

S1_norm = zeros(MC, Nr, Nphi, hsnew);


for i = 1:Nr
    si = squeeze(s(c, i, :));
    si2 = reshape(si, [hs Nphi]).';
    z = zeros(Nphi, hsnew-hs);
    si3 = [si2 z];
    si4 = (reshape(si3.', [1 Nphi*hsnew]));
    
    [S1, F1, Ti1, P1] = spectrogram(si4, hsnew, 0, hsnew, hsnew*sec*1/60); 
    S1_norm(c, i, :, :) = 1./sqrt(hsnew) .* fftshift(S1',2);
%     S1_norm(i, :, :) = S1f(i, :, :)./max(max(squeeze(S1f(i, :, :))));
    S1_norm_db = 20*log10(abs(squeeze(S1_norm(c, i, :, :))));

%     if (n == 1) || (n == length(Vmean))
% 
%     txt = ['SNR = ', num2str(db(SNR(i))/2), ' dB , \Omega = ', num2str(Omega_rpm(m)), ' [rpm]']; %, 'Vmean = ', num2str(Vmean(n))];
% 
%     xl = 'Doppler Velocity [m.s^{-1}]';
%     yl = 'Azimuthal Angle \phi [{\circ}]';
%     zl = 'Power [dB]';
% 
%     surplot(vel_axis_hs, phi_axis*180/pi, S1_norm_db, xl, yl, zl, txt);
%     
%     end
end


%% Calculation of mean velocity and spectrum width


diff_v = diff(vel_axis_hs); dv = diff_v(1);
% PT = zeros(Nr, length(phi_axis));
% mu_re = zeros(Nr, length(phi_axis));
% sigma_re = zeros(1, length(phi_axis));

for i = 1:Nr
    for k = 1:Nphi
        PT_i = squeeze(abs(S1_norm(c, i, k, :))).^2;
        PT(c, i, k, m, n) = sum(PT_i .* dv);
        mu_re_i = vel_axis_hs.' .* squeeze(abs(S1_norm(c, i, k, :))).^2;
        mu_re(c, i, k, m, n) = sum(mu_re_i .* dv)./PT(c, i, k, m, n);
        sigma_re_i = (vel_axis_hs.' - mu_re(c, i, k, m, n)).^2 .* squeeze(abs(S1_norm(c, i, k, :))).^2 .* dv;
        sigma_re(c, i, k, m, n) = sqrt(sum(sigma_re_i)./PT(c, i, k, m, n));
    end
end
end
end

Mu_re_avg = mean(mu_re, 1);
Sigma_re_avg = mean(sigma_re, 1);
Mu_re_std = std(mu_re, 1);
Sigma_re_std = std(sigma_re, 1);



[r_, phi_axis_] = meshgrid(r, phi_axis);

x = r_ .* 1e-3 .* cos(phi_axis_) .* cos(theta); y = r_ .* 1e-3 .* sin(phi_axis_) .* cos(theta);


if (Nr == 1) && (length(Vmean) == 1)
    txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB'];
    dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
    yl = '\mu [m.s^{-1}]';
	xl = 'Azimuthal Angle \phi [{\circ}]';
%     color = 'k';
    marker = markers(m);
    f = figure(103); hold on;
    f.Position = [10 10 1000 1000];
    if Omega_rpm(m) < 20
        
        color = 'k';
    else
        
        color = 'k';
    end
    plott(phi_axis.*180/pi, squeeze(Mu_re_avg(1, :, m, n)), xl, yl, txt, 2, dtext, color, marker)
% else
%     txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB , \Omega , ', num2str(Omega_rpm(m)), ' [rpm]'];
% 
%     xl = 'x [km]';
%     yl = 'y [km]';
%     zl = 'V_{r} Mean [m.s^{-1}]';
%     
%     
%     surplot(x, y, mu_re.', xl, yl, zl, txt); 
end

if (Nr == 1) && (length(Vmean) == 1)
    txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
    dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
    
	xl = 'Azimuthal Angle \phi [{\circ}]';
    
    marker = markers(m);
    f = figure(102); hold on;
    f.Position = [10 10 1000 1000];
     if Omega_rpm(m) < 20
        yyaxis left;
        color = 'k';
        yl =  ['\sigma [m.s^{-1}] for \Omega < 20 [rpm]'];
    else
        yyaxis right;
        color = 'k';
        yl = '\sigma [m.s^{-1}]';
    end
    plott(phi_axis.*180/pi, squeeze(Sigma_re_avg(1, :, m, n)), xl, yl, txt, 2, dtext, color, marker)
% else
%     txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
%     xl = 'x [km]';
%     yl = 'y [km]';
%     zl = 'V_{r} width [m.s^{-1}]';
%     
%     
%     surplot(x, y, sigma_re.', xl, yl, zl, txt);
end

%% Plot with mean Doppler
% for m = 1:length(Omega_rpm)
if length(Vmean) > 1
%Mean of Means for MC rotations with the same ground truth
txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
xl = '\mu True [m.sec^{-1}]';
    
marker = markers(m);
f = figure(104); hold on;   f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['Mean of \mu retrieved [m.s^{-1}], for ', num2str(MC), ' rotations'];
plott(Vmean, squeeze(Mu_re_avg(1, 1, 1, m, :)), xl, yl, txt, 2, dtext, color, marker)

%Mean of Widths for MC rotations with the same ground truth
txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
xl = '\mu True [m.sec^{-1}]';

f = figure(105); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['Mean of \sigma retrieved [m.s^{-1}], for ', num2str(MC), ' rotations'];
plott(Vmean, squeeze(Sigma_re_avg(1, 1, 1, m, :)), xl, yl, txt, 2, dtext, color, marker)

% std of Mean for MC rotations

txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
xl = '\mu True [m.sec^{-1}]';

f = figure(106); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['Std of \mu retrieved [m.s^{-1}], for ', num2str(MC), ' rotations'];
plott(Vmean, squeeze(Mu_re_std(1, 1, 1, m, :)), xl, yl, txt, 2, dtext, color, marker)

% std of Sigma for MC rotations

txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
xl = '\mu True [m.sec^{-1}]';

f = figure(107); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['Std of \sigma retrieved [m.s^{-1}], for ', num2str(MC), ' rotations'];
plott(Vmean, squeeze(Sigma_re_std(1, 1, 1, m, :)), xl, yl, txt, 2, dtext, color, marker)
end
% end
%%
end


%% VAD Technique
Nf = 2;
for n = 1:length(Vspread)
    for i = 1:length(Omega_rpm)
        mu_ret = Mu_re_avg(1, 1, :, i, n);
        [a(i, n, :), b(i, n, :)] = FourierCoeff(mu_ret, Nf, phi_axis);
    end
end