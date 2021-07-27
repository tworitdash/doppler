clear;
%% Ground Truth 
[x, y, z] = meshgrid(linspace(-15e3, 15e3, 100), linspace(-15e3, 15e3, 100), linspace(0, 500, 10));

w_ = 5;

% velocity

H = 10 .* (1 - exp(-z/100));

beta = (40 + 50 .* (1 - exp(-z/100))) .* pi/180;

u = H  .* cos(beta);
v = H .* sin(beta);
w = w_ .* ones(size(x));


%% Radar forward model

r_ = sqrt(x.^2 + y.^2);
phi_ = unwrap(angle(x + 1j .* y));

theta_ = unwrap(angle(x + 1j .* z));

vr = cos(theta_) .* cos(phi_) .* u + cos(theta_) .* sin(phi_) .* v + sin(theta_) .* w;

% figure; surface(squeeze(x(:, :, end)) .* 1e-3, squeeze(y(:, :, end)) .* 1e-3, vr(:, :, end)); shading flat; colormap('jet'); colorbar;

% figure; surface(squeeze(x(:, :, 1)) .* 1e-3, squeeze(y(:, :, 1)) .* 1e-3, squeeze(u(:, :, 1))); shading flat; colormap('jet'); colorbar;

figure; 
h = pcolor(squeeze(x(:, :, end)) .* 1e-3, squeeze(y(:, :, end)) .* 1e-3, vr(:, :, end));
set(h,'ZData',-1+zeros(size(squeeze(vr(:, :, end)))));
shading interp;
colormap('jet');
hold on;
% surface(x_, y_, (vr)); shading flat; colormap('jet'); colorbar;
% figure; surface(r_ .* cos(phi_), r_ .* sin(phi_), (vr.')); shading flat; colormap('jet'); colorbar;
e = 1;
quiver(squeeze(x(1:e:end, 1:e:end, end)), squeeze(y(1:e:end, 1:e:end, end)), squeeze(u(1:e:end, 1:e:end, end)), squeeze(v(1:e:end, 1:e:end, end)), 2, 'w');

figure; plot(squeeze(z(1, 1, :)), squeeze(beta(1, 1, :)) .* 180/pi); grid on;
figure; plot(squeeze(z(1, 1, :)), squeeze(H(1, 1, :))); grid on;