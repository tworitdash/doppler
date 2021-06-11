clear; 
close all;


lambda = 0.03;
M = 360;
hs = 128;

N = hs * M;

mu = 4;


PRT = 1e-3; 

ths = 0:PRT:(hs-1)*PRT;

t = repmat(ths, [1, M]);

t1 = 0:PRT:(N - 1)*PRT;



s = exp(1j * (4 * pi * mu / lambda) .* t);

ph_ = (4 * pi * mu / lambda) .* t1;
s_ = exp(1j .* ph_);



s_f = fftshift(fft(s_));

v_amb = 7.5;

v_axis = linspace(-N/2, N/2-1, N)/N .* 2 * v_amb;

figure; plot(v_axis, db(abs(s_f))); title('Original spectrum')


% th = pi/2;

BW = 1*pi/180;
p0 = 0*pi/180;
p1 = 360*pi/180;

th = linspace(p0, p1, N);

phi = linspace(th(1), th(end), M);

% s_man = s .* exp(1j .* 4*pi/lambda * mu * (cos(th) - 1) .* t);
s_man = abs(s) .* exp(1j .* (angle(s)) + 1j * 4*pi/lambda * mu * (cos(th) - 1) .* t);


s_man_ = abs(s_) .* exp(1j .* (angle(s_) .* cos(th))); 


% + 1j * 4*pi/lambda * mu * (cos(th) - 1) .* t1);% .* exp(1j .* 4*pi/lambda * mu * (cos(th) - 1) .* t1);
% s_man = abs(s_man_) .* exp(1j .* angle(s_man_) .* t./t1);
% s_man_ = s_ .* exp(1j .* 4*pi/lambda * mu * (cos(th) - 1) .* t1);



for i = 1:length(phi)
%     th_(i, :) = th((i - 1)*hs+1:hs*i);
%     t_ = 0:PRT:(hs-1)*PRT;
    s_man_i(i, :) = s_man((i - 1)*hs+1:hs*i);
    s_man_f(i, :) = fftshift(fft(s_man_i(i, :)));

%     figure; plot(v_axis, db(abs(s_man_f))); title('Manipulated spec')
end

vel_axis_hs = linspace(-hs/2, hs/2-1, hs)/hs .* 2 * v_amb;


figure; imagesc(vel_axis_hs, phi*180/pi, db(abs(s_man_f))); title('Manipulated spec'); colormap('jet');


% figure; plot(vel_axis_hs, db(abs(s_man_f(end, :)))); title('Manipulated spec');
% 
% figure; plot(th_(end, :)); title('theta');

s_man_f_full = fftshift(fft(s_man_));

s_man_comp = abs(s_man_) .* exp(1j .* (angle(s_man_)) ./ cos(th));

s_man_comp_f = fftshift(fft(s_man_comp));


figure(1); hold on; plot(v_axis, db(abs(s_man_f_full))); 
hold on; plot(v_axis, db(abs(s_man_comp_f))); 
legend({'HD spectrum', 'Manipulated with cosine Spectrum', 'Compensated from Manipulated spectrum'});  grid on;


%% STFT


