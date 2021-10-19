clear;
close all;

%% Random Fields in 2D


dx = 10;
dy = 10;
x_ = -5e3:dx:5e3;
y_ = -5e3:dy:5e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 


u = 5 .* ones(size(x)); v = ones(size(y));

x0 = [-4e3 -1e3 0 1e3 4e3];
y0 = [-4e3 -1e3 0 1e3 4e3];

sigma_x0 = [5e2 1e3 5e2 1e3 2e2];
sigma_y0 = [5e2 2e3 3e2 2e2 2e3];

R0 = zeros(size(x));

for i = 1:length(x0)
    for k = 1:length(y0)
        R0 = R0 + exp(-((x - x0(i)).^2/(2.*sigma_x0(i).^2) + (y - y0(k)).^2/(2.*sigma_y0(k).^2)));
    end
end

figure; surface(x, y, abs(R0)); shading interp; colorbar; 