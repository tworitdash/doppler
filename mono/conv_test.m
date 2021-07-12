N = 60000; 
PRT = 1e-3;
t = 0:PRT:(N -1)*PRT;
mu = 5;

phi_axis = linspace(0, 2*pi, N);

s = exp(1j .* 4*pi/lambda * mu * sin(phi_axis) ./ (1 * 2*pi/60));
% s = exp(1j .* 4*pi/lambda * 5 * t);
v_amb = 7.5; 

v = linspace(-v_amb, v_amb, N);

s_f = fftshift(fft(s));

figure; plot(v, abs(s_f));

sigma = 0.2;

g_f = 1/(sqrt(2 * pi * sigma^2)) .* exp(-(v - mu).^2./(2 .* sigma^2));

figure; plot(v, abs(g_f))

sg_conv = s_f .* g_f;

figure; plot(v, abs(sg_conv));