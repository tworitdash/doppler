lambda = 0.03;




mu = 4; 
sigma = 0.1;

v_amb = 7.5;

PRT = 1e-3;


BW = 1.8*pi/180;
M = round(2*pi/BW);

Omega = 1*2*pi/60;

hs = round(BW/Omega/PRT);

N = hs .* M;


 th = linspace(0, 2*pi, N);
vel_axis = linspace(-v_amb, v_amb, N);

t = 0:PRT:(N - 1)*PRT;

s = 1/sqrt(2*pi*sigma^2) .* exp(-t .* ((sigma .* 2 /lambda)^2 .* 4 .* pi^2 .* t - 1j .* 2 .* pi .* (mu .* 2./lambda)));

S = fftshift(fft(s));

figure; plot(vel_axis, db(abs(S)));