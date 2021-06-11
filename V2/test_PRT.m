clear; 
close all;


lambda = 0.03;
M = 61;
hs = 128;

N = hs * M;

mu = 4;


PRT = 1e-3; 

ths = 2*hs*PRT:PRT:(3*hs-1)*PRT;

t = repmat(ths, [1, M]);

t1 = 0:PRT:(N - 1)*PRT;

s = exp(1j .* 4 * pi ./ lambda .* mu .* t);

s_ = exp(1j .* 4 * pi ./ lambda .* mu .* t1);

s_f = fftshift(fft(s_));

v_amb = 7.5;

v_axis = linspace(-N/2, N/2-1, N)/N .* 2 * v_amb;

figure; plot(v_axis, db(abs(s_f))); title('Original spectrum')
