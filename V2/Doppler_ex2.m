clear; 
close all;


lambda = 0.03;
M = 4;
hs = 16;

N = hs * M;

mu = 4;

v_amb = 7.5;

PRT = 1e-3; 


t1 = 0:PRT:(N - 1)*PRT;


% ths = 0:PRT:(hs-1)*PRT;
% 
% t = repmat(ths, [1, M]);
% % 

ph_ = (4 * pi * mu / lambda) .* t1;
s_ = exp(1j .* ph_);

p0 = 0*pi/180;
p1 = 360*pi/180;

th = linspace(p0, p1, N);
% th = pi/3;
% th = 0;

phi = linspace(th(1), th(end), M);

% s_man = abs(s_) .* exp(1j * (angle(s_))  + 1j * (4*pi/lambda * mu * (cos(th) - 1) .* t1));

% s_man = s_ .* exp(1j * (4*pi/lambda * mu * (cos(th) - 1) .* t1));

s_man = abs(s_) .* exp(1j .* unwrap(angle(s_)) .* sin(th) ./ th);




for i = 1:length(phi)
    th_ = th((i - 1)*hs+1:hs*i);
    s_man_i(i, :) = abs(s_((i - 1)*hs+1:hs*i)) .* exp(1j * unwrap(angle(s_((i - 1)*hs+1:hs*i))) .* cos(th_));
    
%     s_man_i(i, :) = s_man((i - 1)*hs+1:hs*i) ;
    s_man_f(i, :) = fftshift(fft(s_man_i(i, :)));
    
end

s_man_i_reshape = reshape(s_man_i, 1, N);

vel_axis_hs = linspace(-hs/2, hs/2-1, hs)/hs .* 2 * v_amb;


figure; imagesc(vel_axis_hs, phi*180/pi, db(abs(s_man_f))); title('Manipulated spec'); colormap('jet');


%% STFT



[S1, F1, Ti1, P1] = spectrogram(s_man, hs, 0, hs, N/PRT);


S1 = S1./max(S1(:));
S1_norm_db = 20*log(abs(fftshift(S1',2)));

F1_doppler = F1 - N/2;

V_doppler = lambda .* F1_doppler ./ 2;

figure;
imagesc(vel_axis_hs, phi * 180/pi, S1_norm_db);
colormap('jet');
colorbar;

xlabel('Doppler velocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Azimuth [Deg]', 'FontSize', 12, 'FontWeight', 'bold');
%title([' Micro Doppler N_{fft} = ', num2str(nfft), ', overlap = ', num2str(hop),...
%    ', fs = ', num2str(fs*1e-6), ' MHz, ', 'WINDOW = ', num2str(wlen), ' er = 5'], 'FontSize', 10, 'FontWeight', 'bold');
% title([' Micro Doppler N_{fft} = ', num2str(nfft), ', overlap = ', num2str(hop),...




%% Phase change


figure; plot(diff(angle(s_man_i_reshape))); hold on; plot(diff(angle(s_man)));