clear;
% close all;
%% Generating time domain data for a monochriomatic wind within a rotating radar

BW_deg = 1;
n_rot = 1;

phi_0_deg = 0;
PRT = 1e-3;
Omega_rpm = 6;

lambda = 0.03;
%% 

BW = BW_deg * pi/180;
Omega = Omega_rpm * 2*pi/60; % rotation speed in rad/s
Td = BW/Omega;
hs = round(Td/PRT);
phi_0 = phi_0_deg * pi/180;

sec = round((n_rot*2*pi)/BW);
N = sec * hs;
t = 0:PRT:(N - 1)*PRT;
% phi_axis = phi_0:BW:phi_0+Omega*t(end);
phi_axis = linspace(phi_0, Omega*t(end), sec);
% phi_axis = zeros(size(phi_axis));


%% variables for signal model

beta_wind_deg = 0;

beta_wind = beta_wind_deg .* pi/180;

mu = normrnd(3, 0.1, [1 N]);
n_sig = 0;
% mu = 3;
[s] = TD_generator(mu, lambda, beta_wind, phi_0, Omega, t, n_sig);

% figure; plot(t, real(s)); hold on; plot(t, imag(s));

v_amb = lambda/(4 * PRT);

% if mod(N, 2) == 0
%     vel_axis = linspace(-N/2, N/2-1, N)./N .* 2 .* v_amb;
% else
%     vel_axis = linspace(-v_amb, v_amb, N);
% end

% figure; plot(vel_axis, db(abs(fftshift(fft(s))))); grid on;


%% Pulse-Pair method


if mod(hs, 2) == 0
    vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
else
    vel_axis_hs = linspace(-v_amb, v_amb, hs);
end
s_seg = reshape(s, hs, sec);


k = 1:1:round(hs/2);
% k = 1;
R = zeros(sec, length(k));

for l = 1:length(k)
    for i = 1:sec
        R(i, l) = 1./hs .* [s_seg(1:end-k(l), i)'] * [s_seg(1+k(l):end, i)];
        fn(i, l) = 1./(2 .* pi .* PRT .* k(l)) .* angle(R(i, l));
    end


end

figure(101); plot(phi_axis .* 180/pi, fn .* lambda/2, 'LineWidth', 2); grid on;

for i = 1:sec

    f_vn(i) = (1/(2 * pi * PRT)) .* angle(sum(abs(R(i, :)) .* exp(1j .* 2 * pi .* fn(i, :) .* k .* PRT)));
    
    f_vun(i) = (1/(2 * pi * PRT)) .* angle(sum(abs(R(i, :)).^(1./k) .* exp(1j .* 2 * pi .* fn(i, :) .* k .* PRT)));
    
    f_av(i) = 1/length(k) .* sum(fn(i, :));
end

figure; plot(phi_axis .* 180/pi, f_vn .* lambda/2, 'LineWidth', 2); grid on;
hold on; plot(phi_axis .* 180/pi, f_vun .* lambda/2, 'LineWidth', 2); 
hold on; plot(phi_axis .* 180/pi, f_av .* lambda/2, 'LineWidth', 2);
