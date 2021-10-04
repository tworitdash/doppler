clear;
close all;


dx = 1e-3;
dy = 1e-3;

x_ = eps:dx:10;
y_ = eps:dy:3;

R0 = rand(length(x_), (length(y_))) + 1j .* rand(length(x_), (length(y_))); 


figure; surface(x_, y_, abs(R0).'); shading flat;
F = fft2(R0);

fx = 0.005;
fy = 0.002; 

nx1 = 1 + floor(0.9999.*fx*length(x_));
nx2 = length(x_) - nx1 + 1; 

ny1 = 1 + floor(0.9999.*fy*length(y_));
ny2 = length(y_) - ny1 + 1;


F(nx1:nx2, :) = 0;
F(:, ny1:ny2) = 0;

R1 = ifft2(F);
% R1            = imgaussfilt(real(F), 0.001) + 1j .* imgaussfilt(imag(F), 0.0011);

figure; surface(x_, y_, abs(R1).'); shading flat; colorbar; colormap('jet')