clear;
close all;
c0 = 3e8;

PRT = 0.01;


N = 3;


t_ = 0:PRT:(N - 1)*PRT;

Rmax = 25e3;

dR = 2000;

x_ = -Rmax:dR:Rmax;
y_ = -Rmax:dR:Rmax; 

[x, y, t] = meshgrid(x_, y_, t_);

D = [1 0; 0 1];

w_ = [1 0].';

alpha = 0.33;

sigma_s = 1e-3;

gamma = 333;
w = [0 0].';
for i = 1:size(x, 1)
    for k = 1:size(x, 2)
%        q(i, k, :) = sigma_s .* exp(-sqrt((x(1, k, 1) - x(1, i, 1)).^2 + (y(k, 1, 1) - y(i, 1, 1)).^2)./gamma) .* randn(1, size(x, 3));
%          q(i, k, :) =  normrnd(0, sqrt(sigma_s), [1 size(x, 3)]);
         w = w + PRT .* w_;
         q(i, k, :) = sqrt(sigma_s) .* rand([1 size(x, 3)]);
         H(i, k) = alpha .* exp(-([x(1, k, 1) y(k, 1, 1)] - [x(1, i, 1) y(i, 1, 1)] - w.') * inv(D) * ([x(1, k, 1) y(k, 1, 1)] - [x(1, i, 1) y(i, 1, 1)] - w.').');
    end
end

[R, C] = meshgrid(x_, y_);

u = gauss2d(R, C, 5000, [x_(14) y_(14)]);


figure; surface(x_, y_, u); shading flat; colormap('jet');

u_t = zeros(size(x));

u_t(:, :, 1) = u;

for l = 2:size(x, 3)

    u_t(:, :, l) = H * u_t(:, :, l - 1); % + q(:, :, l);
    figure; surface(x_, y_, squeeze(u_t(:, :, l))); shading flat; colormap('jet');
%     figure; surface(x_, y_, squeeze(q(:, :, l))); shading flat; colormap('jet');
end



