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
Omega_rpm = 60;
% Omega_rpm = linspace(1, 60, 4);
% Omega_rpm = 60;

lambda = 0.03;



BW = BW_deg * pi/180;
Omega = Omega_rpm * 2*pi/60; % rotation speed in rad/s
for m = 1:length(Omega_rpm)
    
clear s; 


% T = (2*pi)/Omega(m);

% N = ceil(T/PRT);

% t = 0:PRT:(N - 1)*PRT;

sec = ceil(2*pi/BW);

ph = 0:BW:(sec)*BW;

T = ph(end)/Omega(m);


N = ceil(T/PRT);

hs = ceil(N/sec);

N = hs*sec;
t = 0:PRT:(N - 1)*PRT;


% Td = BW/Omega(m);
% hs = round(Td/PRT);

hs_(m) = hs;

phi_0 = phi_0_deg * pi/180;
% phi_axis_1 = phi_0:BW:phi_0+Omega(m)*t(end);
phi_axis = linspace(phi_0, phi_0+Omega(m)*t(end), sec);

% phi_axis = mean([phi_axis_1(1:end-1); phi_axis_1(2:end)]);
Nphi = length(phi_axis);
% phi_axis = zeros(size(phi_axis));

%% variables for signal model

% beta_wind_deg = [0 .* ones(1, round(N/2)) 90 .* ones(1, N - round(N/2))];
% beta_wind_deg  = 3 .* Omega .* t .* 180/pi;
beta_wind_deg = 0;
theta = 0;

beta_wind = beta_wind_deg .* pi/180;

v_amb = lambda/(4 * PRT);

% mu_mean = linspace(-v_amb, v_amb, 100);

% Vmean = normrnd(5, 1, [1 Nr]);

Vmean = 3;

% u_mean = normrnd(Vmean .* cos(beta_wind), 0.1, [1 Nr]);
% v_mean = normrnd(Vmean .* sin(beta_wind), 0.1, [1 Nr]);

% u_phi = zeros(Nphi, Nr); v_phi = zeros(Nphi, Nr); 
u_phi = 0; v_phi = 0;

% pd = makedist('Weibull');
s = zeros(Nr, N);
for i = 1:Nr
   
%     mu = normrnd(mu_mean(i), 0.1, [1 N]);
%     u = normrnd(u_mean(i), 0.2, [1 N]);
%     v = normrnd(v_mean(i), 0, [1 N]);
    
    V = normrnd(Vmean(i), 0.1, [1 N]);
%     V = random('burr', Vmean(i), 2, 5, 1, N);
    u = V .* cos(beta_wind);
    v = V .* sin(beta_wind);
    
%     V = sqrt(u.^2 + v.^2);
    
    mu = V .* cos((beta_wind - Omega(m) .* t - phi_0)) .* cos(theta); % +  v .* sin(beta_wind - Omega .* t - phi_0) .* cos(theta);
    
    p = Omega(m) .* t *180/pi;
    figure(103); hold on; plot(p, mu, 'DisplayName', [num2str(Omega_rpm(m))]); legend;
%     u_reshape = reshape(u, [Nphi hs]);
%     u_phi(:, i) = mean(u_reshape, 2);
%     
%     v_reshape = reshape(v, [Nphi hs]);
%     v_phi(:, i) = mean(v_reshape, 2);
    
%     mu = u.' .* cos(beta_wind - Omega .* t - phi_0) +  v.' .* sin(beta_wind - Omega .* t - phi_0);

    % figure; quiver(cos(Omega .* t), sin(Omega .* t), u, v);

    % mu = normrnd(4, 1, [1 N]);
    n_sig = sqrt(0.001); 
    % n_sig = 0;

    % mu = 3;
    [s(i, :), SNR(i)] = TD_generator(mu, lambda, beta_wind, phi_0, Omega(m), t, n_sig);
end


%% Retrieval

% if (mod(hs, 2) == 0) && (hs < 16)
%     hsnew = hs + 1;
% %     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
%      vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
% elseif (mod(hs, 2) ~= 0) && (hs < 16)
%     
%     if (2^(nextpow2(hs)) - hs) == 3
%         hsnew = hs+2;
%         vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
%     else
%         hsnew = hs;
%         vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
%     end
% elseif (mod(hs, 2) == 0) && (hs > 16)
%      hsnew = hs;
% %     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
%      vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
% elseif (mod(hs, 2) ~= 0) && (hs > 16)
%     hsnew = hs;
%     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
% end

% ================================================================================
% if (mod(hs, 2) == 0) && (hs > 16)
%     hsnew = hs;
%     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
% elseif (hs < 16) && (mod(hs, 2) == 0)
% %     hsnew = 2^(nextpow2(hs)) + 1;
%     hsnew = hs;
% %     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
%     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
% else
%     hsnew = hs;
%     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
% end

% ================================================================================
% 
% if (mod(hs, 2) == 0) && (hs > 16)
%     hsnew = hs;
%     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
% elseif (hs < 16) 
% %     hsnew = 2^(nextpow2(hs)) + 1;
%     hsnew = hs;
% %     hsnew = 64;
%     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
% %     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
% else
%     hsnew = hs;
%     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
% end

if (mod(hs, 2) == 0)
    hsnew = hs;
    vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
%     vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
else
    hsnew = hs;
    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
end




S1_norm = zeros(Nr, Nphi, hsnew);

for i = 1:Nr
    si = s(i, :);
    si2 = reshape(si, [hs Nphi]).';
    z = zeros(Nphi, hsnew-hs);
    si3 = [si2 z];
    si4 = (reshape(si3.', [1 Nphi*hsnew]));
    
    [S1, F1, Ti1, P1] = spectrogram(si4, hsnew, 0, hsnew, hsnew*sec*1/60); 
    S1_norm(i, :, :) = 1./sqrt(hsnew) .* fftshift(S1',2);
%     S1_norm(i, :, :) = S1f(i, :, :)./max(max(squeeze(S1f(i, :, :))));
    S1_norm_db = 20*log10(abs(squeeze(S1_norm(i, :, :))));


%     txt = ['SNR = ', num2str(db(SNR(i))/2), ' dB , \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
% 
%     xl = 'Doppler Velocity [m.s^{-1}]';
%     yl = 'Azimuthal Angle \phi [{\circ}]';
%     zl = 'Power [dB]';
% 
%     surplot(vel_axis_hs, phi_axis*180/pi, S1_norm_db, xl, yl, zl, txt);
    
end

%% Calculation of mean velocity and spectrum width


diff_v = diff(vel_axis_hs); dv = diff_v(1);
PT = zeros(1, length(phi_axis));
mu_re = zeros(1, length(phi_axis));
sigma_re = zeros(1, length(phi_axis));

for i = 1:Nr
    for k = 1:Nphi
        PT_i = squeeze(abs(S1_norm(i, k, :))).^2;
        PT(i, k, m) = sum(PT_i .* dv);
        mu_re_i = vel_axis_hs.' .* squeeze(abs(S1_norm(i, k, :))).^2;
        mu_re(i, k, m) = sum(mu_re_i .* dv)./PT(i, k, m);
        sigma_re_i = (vel_axis_hs.' - mu_re(i, k, m)).^2 .* squeeze(abs(S1_norm(i, k, :))).^2 .* dv;
        sigma_re(i, k, m) = sqrt(sum(sigma_re_i)./PT(i, k, m));
    end
end

[r_, phi_axis_] = meshgrid(r, phi_axis);

x = r_ .* 1e-3 .* cos(phi_axis_) .* cos(theta); y = r_ .* 1e-3 .* sin(phi_axis_) .* cos(theta);


if Nr == 1
    txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB , \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
    dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
    yl = 'Mean Doppler Velocity [m.s^{-1}]';
	xl = 'Azimuthal Angle \phi [{\circ}]';
    color = 'k';
    marker = markers(m);
    figure(103); hold on;
    plott(phi_axis.*180/pi, squeeze(mu_re(1, :, m)), xl, yl, txt, 2, dtext, color, marker)
else
    txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB , \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];

    xl = 'x [km]';
    yl = 'y [km]';
    zl = 'V_{r} Mean [m.s^{-1}]';
    
    
    surplot(x, y, mu_re.', xl, yl, zl, txt); 
end

if Nr == 1
    txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
    dtext = [' \Omega = ', num2str(Omega_rpm(m)), ' [rpm]'];
    yl = 'Doppler spectrum width [m.s^{-1}]';
	xl = 'Azimuthal Angle \phi [{\circ}]';
    color = 'k';
    marker = markers(m);
    figure(102); hold on;
    plott(phi_axis.*180/pi, squeeze(sigma_re(1, :, m)), xl, yl, txt, 2, dtext, color, marker)
else
    txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB '];
    xl = 'x [km]';
    yl = 'y [km]';
    zl = 'V_{r} width [m.s^{-1}]';
    
    
    surplot(x, y, sigma_re.', xl, yl, zl, txt);
end
end

%% Calculation of wind fields [VAD Technique]
% Vr = @(a0, a1, b1, a2, b2) a0 + a1 .* cos(phi_axis) + a2 .* cos(2.*phi_axis) + ...
%     b1 .* cos(phi_axis) + b2 .* cos(2.*phi_axis);
Nf = 2;
for i = 1:Nr

%     param = fit(phi_axis.', mu_re(i, :).', 'fourier2');
%     a0(i) = param.a0;
%     a1(i) = param.a1;
%     b1(i) = param.b1;
%     a2(i) = param.a2;
%     b2(i) = param.b2;
%     w(i) = param.w;

    [a(i, :), b(i, :)] = FourierCoeff(mu_re(i, :), Nf, phi_axis);
%     figure; plot(param, phi_axis, mu_re(i, :));
%     Vfs_re(:, i) = a0(i) + a1(i) .* cos(w(i) .* phi_axis) +  b1(i) .* sin(w(i) .* phi_axis)...
%         + a2(i) .* cos(2 .* w(i) .* phi_axis) +  b2(i) .* sin(2.*w(i) .* phi_axis);
    
end
Vfs_re = zeros(Nphi, Nr);
for i = 1:Nr
    Vfs_re(:, i) =  a(i, 1) .* ones(1, Nphi);
    for k = 2:Nf+1
        Vfs_re(:, i) = Vfs_re(:, i) + (a(i, k) .* cos([k - 1].*phi_axis)).' + (b(i, k) .* sin([k - 1].*phi_axis)).';
    end
end

a1 = a(:, 2); b1 = b(:, 2);
a2 = a(:, 3); b2 = b(:, 3);

V_vad = sqrt(a1.^2 + b1.^2)./(cos(theta));
V_vaddeform = 2*sqrt(a2.^2 + b2.^2)./(r.' .* cos(theta));

phi_vad = zeros(1, Nr);
phi_vaddeform = zeros(1, Nr);

for i = 1:Nr
    if b1(i) < 0
        phi_vad(i) = pi/2 - atan2(a1(i),b1(i));
    else
        phi_vad(i) = pi/2 - atan2(a1(i),b1(i));
    end
    if b2(i) < 0
        phi_vaddeform(i) = pi/4 - 1/2 * atan2(a2(i),b2(i));
    else
        phi_vaddeform(i) = pi/4 - 1/2 * atan2(a2(i),b2(i));
    end
    
end

V_vadph = repmat(V_vad, 1, Nphi);
V_vadph = reshape(V_vadph, [Nr Nphi]);
phi_vadph = repmat(phi_vad, 1, Nphi);
phi_vadph = reshape(phi_vadph, [Nr Nphi]);
V_rvadph = V_vadph.' .* cos(phi_vadph.' - phi_axis_ - phi_0);

V_vaddeformph = repmat(V_vaddeform, 1, Nphi);
V_vaddeformph = reshape(V_vaddeformph, [Nr Nphi]);
phi_vaddeformph = repmat(phi_vaddeform, 1, Nphi);
phi_vaddeformph = reshape(phi_vaddeformph, [Nr Nphi]);



txt = ['Wind fields from first harmonic'];
% 
xl = 'x [km]';
yl = 'y [km]';
ex = 10; ey = 1;

figure; quiverr(x, y, u_phi, v_phi, xl, yl, ex, ey, txt);
hold on; quiverr(x, y, V_vadph.'.*cos(phi_vadph).' ,  V_vadph.'.*sin(phi_vadph).', xl, yl, ex, ey, txt);
% figure; plot(param, phi_axis, mu_re);

% txt = ['Wind fields from second harmonic'];
% % 
% xl = 'x [km]';
% yl = 'y [km]';
% ex = 10; ey = 1;
% 
% figure; quiverr(x, y, u_phi, v_phi, xl, yl, ex, ey, txt);
% hold on; 
figure;
quiverr(x, y, V_vaddeformph.'.*cos(phi_vaddeformph).' , ...
    V_vaddeformph.'.*sin(phi_vaddeformph).', xl, yl, ex, ey, txt);

% txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB , \Omega = ', num2str(Omega_rpm), ' [rpm]'];
txt = ['V_{r} mean VAD [m.s^{-1}] only from one harmonic'];
xl = 'x [km]';
yl = 'y [km]';
zl = 'V_{r} mean VAD [m.s^{-1}] only from one harmonic';


surplot(x, y, V_rvadph, xl, yl, zl, txt);

% txt = ['SNR = ', num2str(db(mean(SNR))/2), ' dB , \Omega = ', num2str(Omega_rpm), ' [rpm]'];
txt = ['V_{r} mean VAD [m.s^{-1}] by fitting the Fourier Series'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'V_{r} mean VAD [m.s^{-1}] by fitting the Fourier Series';


surplot(x, y, Vfs_re, xl, yl, zl, txt);
%     
%%

% Theta = 0:pi/180:2*pi;
% Mu_mean = zeros(1, length(Theta));
% diffphi = diff(phi_axis); dphi = diffphi(1);
% 
% % for k = 1:length(Theta)
% %     Mu_mean(k) = sum(mu_re .* exp(1j .* 2 * pi * cos(phi_axis - Theta(k))) .* dphi);
% % end
% 
% Mu_mean = 1./sqrt(length(phi_axis)) .* (fft(mu_re));
% 
% yl = 'Spatial angular function [m.s^{-1}]';
% xl = '\Theta [^{\circ}]';
% 
% plott(phi_axis*180/pi, abs(Mu_mean).^2, xl, yl, txt, 4);


%% Errors

phiError = (beta_wind - phi_vad) .* 180/pi;
