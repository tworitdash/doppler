v = 4;
lambda = 0.03;
PRT = 1e-3;

t = eps:PRT:1024*PRT;

dis = v .* t;
% 
for i = 1:length(t)
    if mod(i, 5) == 0
        dis(i) = 0;
    end
end

s = exp(1j .* 4 * pi/lambda * dis);

sfft = 1./sqrt(length(s)) .* fftshift(fft(s));

v_amb = lambda/(4 * PRT);

v_axis = linspace(-v_amb, v_amb, length(t));

figure; plot(v_axis, db(abs(sfft)));

