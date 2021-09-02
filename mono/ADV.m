%% Advection models for clouds
clear; 
close all;

%% Define the space

N = 20;

delx = 250;
dely = 250;

x = linspace(0, delx*N, N);
y = linspace(0, dely*N, N);

[X, Y] = meshgrid(x, y);
Z = rand(size(X));
figure; contour(X, Y, Z);

% Example - syntehsis of a likely Digital Elevation Model
H=2; 

dem = gen_UMF2d(1.8, 0.05, H, size(X, 1)); 
colormap('bone'); surfl(dem); shading('interp')


