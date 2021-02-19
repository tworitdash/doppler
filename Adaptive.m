% Script goal:
% Simulate a signal that contains phase information

% Script options:
% Set the width distribution of the target
% Add some constant velocity to target and check Doppler velocity 
clc
clear
close all

% Radar parameters:
lambda = 0.03; %m
PRT = 1e-3;
hits_scan = 2^20;

delta_v = lambda/(2*hits_scan*PRT);% velocity resolution
v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
vel = linspace(-v_amb,v_amb,hits_scan);% velocity axis

dR = 1;% m % range resolution [m]
Rmax = 1.5e3; %max unabiguous range [m]
Nr = Rmax/dR; % number of range bins
range = (0:Nr-1)*dR; % range bins

% target_pos_0 = [1e3 5e3];% target position in [m]
N = 1; % 

target_pos = abs(randn(1, N)) .* 1e2
target_vel = abs(randn(1, N)) .* 2
        

% For distributed target
tau = .1e3;% 1km target width (to define target distribution) 

% For point target
% tau = .1e3;

% vt = 3;% target velocity [m/s] 
data = [];

pn = 4;
mn = 2;


for i = 1:hits_scan
signal = zeros(1,numel(range));
%1. Update target position
    for k = 1:N
        target_pos_ = target_pos(k) + target_vel(k) * PRT * (i - 1); % + (PRT * i) ** (1) * 1
        target_range = target_pos_ + [0:dR:tau-1];
        ph = 2*pi*2*target_range/lambda;
        target_range_index = round(target_range/dR);
        signal(target_range_index) = 20 .* exp(1j * ph) + randn(1) * 6;
    end
% target_pos = target_pos_0+vt*PRT*(i-1)+1./2 * vt * (PRT * (i - 1)).^2;
% 
% %2. Compute phase w.r.t. to each updated range bin value saved in target_range
% target_range = [0:dR:tau-1]+target_pos;
% ph = 2*pi*2*target_range/lambda; % phase as function of range 2pi*2range/lambda
% 
% %3. Compute indexes that correspond with new range bins of the target
% target_range_index = round(target_range/dR);
% 
% %4. Create the signal with phase information at corresponding range bin value
% 
% signal(target_range_index) = 40 .* exp(1i*ph) + pn .* randn(1);
% %5. Save signals in a matrix
data(i,:) = signal + mn .* randn(1, Nr);
end


a = zeros(hits_scan, Nr);
Theta_model = ones(hits_scan, Nr);
data_model = zeros(hits_scan, Nr);
data_model_var = ones(hits_scan, Nr);
gamma = zeros(hits_scan, Nr);

gamma(1, :) = abs(data(1, :)).^(-2);
P = ones(hits_scan, Nr);
Q = zeros(hits_scan, Nr);

% 
% data_model(1, :) = data(1, :);

tic;

for  k = 2:hits_scan
    %disp(k)
    
    gamma(k, :) = gamma(k - 1, :) ./ (1 + gamma(k - 1, :) .* abs(data_model(k - 1, :)).^2);
    a(k, :) = a(k - 1, :) + gamma(k, :) .* conj(data_model(k -1, :)) .* (data(k, :) - a(k - 1, :) .* data_model(k - 1, :));
    
    
    
    Q(k - 1, :) = sum(abs(data(1:k-1, :)), 1).^2 ./ (mn .* sum(1:k - 1) - 1);
  
    P(k, :) = (Q(k - 1, :) .* (1 - abs(a(k, :)).^2) + P(k - 1, :) .* abs(a(k, :)).^2) ./ ...
        (1 + Q(k - 1, :) .* (1 - abs(a(k, :)).^2) + P(k - 1, :) .* abs(a(k)).^2);
    data_model(k, :) = a(k, :) .* data_model(k - 1, :) + P(k, :) .* (data(k, :) - a(k, :) .* data_model(k - 1, :));
    
    
end
omega = angle(a);

F = pn.^2./(1 - a .* exp(-1j .* omega));
v = F .* lambda / 2;
time_adap = toc;


v_max = 10;
v_abs = abs(v);
v_abs(v_abs>v_max) = v_max;

data_model_abs = abs(data_model);
data_model_max = 70;
data_model_abs(data_model_abs>data_model_max) = data_model_max;

figure; surface(range*1e-3, 1:hits_scan, abs(data)); shading flat; colorbar, colormap('jet'); 

xlabel('Range [km]')
ylabel('Hits scan')
title('Raw data in time domain')

figure; surface(range*1e-3, 1:hits_scan, data_model_abs); shading flat; colorbar, colormap('jet');caxis([0 70]);

xlabel('Range [km]')
ylabel('Hits scan')
title('Data after filtering, Estimate of the compelx data')

figure; surface(range*1e-3, 1:hits_scan, v_abs); shading flat; colorbar, colormap('jet'); caxis([1 10]);

xlabel('Range [km]')
ylabel('Hits scan')
title('velocity of the targets [m/s]')


figure;
imagesc(range*1e-3,[],abs(data)); colorbar;
xlabel('Range [km]')
ylabel('Hits per scan')
title('Raw data 2D plot')

%6. Make Doppler
tic;
data_Doppler = fftshift(fft(data,[],1), 1);
time_fft = toc;
figure;
imagesc(vel,range*1e-3,abs(data_Doppler).'); shading flat; colorbar; colormap('jet');
% surface(vel,range*1e-3,abs(data_Doppler).'); shading flat; colorbar; colormap('jet');
xlabel('Velocity [m/s]')
ylabel('Range [km]')
title('Range Doppler Map')