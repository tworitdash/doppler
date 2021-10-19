clear;

%% velocity ground truth

BW = 1.8*pi/180;
phi_axis = eps:BW:2*pi-BW;
PRT = 1e-3;

beta_wind = eps;
N = 5;
u = 5;
n_rot = 1000;


sec = round((2*pi)/BW);
for i = 1:length(phi_axis)-1
    phi_all((i - 1)*N+1:i*N) = linspace(phi_axis(i), phi_axis(i+1), N);
end

vr = u .* cos(phi_all - beta_wind);


for k = 1:length(phi_all)
    t = [(k - 1) * PRT (k - 1)*PRT+N*sec*(1:n_rot-1)*PRT]; 
    dis(k, :) = vr(k) .* t;
end

%% Radar Doppler signal


lambda = 0.03;
PRT = 1e-3;
v_amb = lambda/(4 * PRT);

N_rot = 4; % number of rotations needed to integrate

v_axis = linspace(-v_amb, v_amb, N_rot*N);

for i = 1:length(phi_axis)-1
    phase_sig(i, :) = exp(-1j .* 4*pi/lambda .* reshape(dis((i - 1)*N+1:i*N, 1:N_rot), [1 N_rot*N]));
end

sfft = 1./sqrt(N_rot*N) .* fftshift(fft(phase_sig, [], 2));


figure; surface(phi_axis(1:end-1)*180/pi, v_axis, db(abs(sfft.'))); shading flat; colormap('jet');

figure; plot(v_axis, db(sfft(end, :)));
