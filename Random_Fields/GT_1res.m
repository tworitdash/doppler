%%
clear;
close all;
markers = load('../mono/markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;

N = 1000;

r0 = 100;
dr = 50;
ph0 = 0*pi/180;
dph = 1.8*pi/180;

r = r0 - dr/2 + dr .* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);

% x0 = 10 + 5 .* rand(1, N);
% % y0 = 10 + 5 .* rand(1, N); 
% y0 = zeros(1, N);

x0 = r .* cos(ph); y0 = r .* sin(ph);


% x0 = linspace(50, 100, N);
% y0 = zeros(1, N);
x(1, :) = x0; y(1, :) = y0;

%% Plot positions of scatterers
% txt = ['Resolution cell'];
% dtext = ['Scaterer position at time t = 0'];
% xl = 'x[m]';
% 
% f = figure(107); hold on; f.Position = [10 10 1000 1000];
% color = 'k';
% yl =  ['y[m]'];
% 
% marker = markers(2);

% plott2(x0, y0, xl, yl, txt, 2, dtext, color, marker)
%% Update with time

x(1, :) = x0;
y(1, :) = y0;

% v = linspace(1, 5, N);
u = normrnd(10, 5, [1 N]);
v = normrnd(0, 0, [1 N]);
Nt = 128;

A0 = ones(1, N);
D0 = 0;

% lambda = 0.03; 

dt = 1e-3;

% z(1) = sum(A0 .* exp(1j .* 4 * pi / lambda * x0));

% z(1) = 0;
% SNR_db = 30; 
% SNR = 10^(SNR_db/10);
 
for i = 2:Nt
    x(i, :) = x(i - 1, :) + u .* dt;
    y(i, :) = y(i - 1, :) + v .* dt;
    
%     figure(101); hold on; plot(x(i, :), y(i, :), '+');
    
    D(i, :) = sqrt(x(i, :).^2 + y(i, :).^2);
%     z(i) = sum(A0 .* exp(1j .* 4 * pi / lambda .* D(i, :)));
%     close all;
%     plott2(x(i, :), y(i, :), xl, yl, txt, 2, dtext, color, marker)
    
end

%% 
figure; 
for i = 1:Nt
     plot(x(i, :), y(i, :), '*');
     pause(0.1);
end

 