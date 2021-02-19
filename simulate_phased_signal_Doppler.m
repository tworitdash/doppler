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
hits_scan = 256;

delta_v = lambda/(2*hits_scan*PRT);% velocity resolution
v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
vel = linspace(-v_amb,v_amb,hits_scan);% velocity axis

dR = 1;% m % range resolution [m]
Rmax = 15e3; %max unabiguous range [m]
Nr = Rmax/dR; % number of range bins
range = (0:Nr-1)*dR; % range bins

target_pos = 5e3;% target position in [m]

% For distributed target
% tau = 1e3;% 1km target width (to define target distribution) 

% For point target
tau = .1e3;

vt = 5.4;% target velocity [m/s] 
data = [];
for i = 1:hits_scan

%1. Update target position
target_pos = 5e3+vt*PRT*(i-1);

%2. Compute phase w.r.t. to each updated range bin value saved in target_range
target_range = [0:dR:tau-1]+target_pos;
ph = 2*pi*2*target_range/lambda; % phase as function of range 2pi*2range/lambda

%3. Compute indexes that correspond with new range bins of the target
target_range_index = round(target_range/dR);

%4. Create the signal with phase information at corresponding range bin value
signal = zeros(1,numel(range));
signal(target_range_index) = exp(1i*ph);

%5. Save signals in a matrix
data(i,:) = signal;
end

figure;
imagesc(range*1e-3,[],abs(data)); colorbar;
xlabel('Range [km]')
ylabel('Hits per scan')
title('Raw data 2D plot')

%6. Make Doppler
data_Doppler = fftshift(fft(data,[],1), 1);

figure;
imagesc(vel,range*1e-3,abs(data_Doppler).'); colorbar;
xlabel('Velocity [m/s]')
ylabel('Range [km]')
title('Range Doppler Map')