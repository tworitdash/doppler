%% 
clear; 
close all;

lambda = 0.03;
%% Radar forward model 

% th = linspace(2.5, 60, 15) * pi/180;
ph = linspace(eps, 360, 360) * pi/180; % 1 degree azimuth resolution
% ph = linspace(eps, 1.88, 1024) * pi/180;

dR = 227;
R_max = 15e3;

r = eps:dR:R_max; % range axis

[R, phi] = meshgrid(r, ph);


P = zeros(size(phi));

P(30:60, :) = 1;

% for i = 30:60
%     P(i, )
% end

figure; surface(R .* cos(phi), R .* sin(phi), P); shading flat; colorbar;