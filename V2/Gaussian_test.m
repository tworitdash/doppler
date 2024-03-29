clear;
close all;

PRT = 1e-3; 

N = 129; 

lambda = 0.03;

v_amb = 1/(PRT);

vel_axis = linspace(-v_amb, v_amb, N);

mu = 4 .* 4 ./ lambda;
sigma = 0.2 .* 4 ./ lambda;

S = 1./(sqrt(2 * pi * sigma^2)) .* exp(-(vel_axis - mu ).^2./(2 * sigma^2));

figure; plot(vel_axis, abs(sqrt(S)));

% figure; plot(vel_axis, real(sqrt(S))); hold on; plot(vel_axis, imag(sqrt(S)));

s = ifftshift(ifft(sqrt(S)));



t = 0:PRT:(N - 1)* PRT;

s_analyt = (2/pi)^(1/4) * sqrt(sigma) .* exp(-t .* (sigma^2 .* t + 1j .* mu));

figure; plot(t, real(s)); hold on; plot(t, imag(s));
figure; plot(t, abs(s));

figure; plot(t, real(s_analyt)); hold on; plot(t, imag(s_analyt));
figure; plot(t, abs(s_analyt));

