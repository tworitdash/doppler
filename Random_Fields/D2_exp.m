%% 2D echo simulator test
clear; close all;
markers = load('../mono/markers.mat');
markers = markers.markers;


x_ = linspace(eps, 30, 200);
y_ = linspace(eps, 30, 200);

[x, y] = meshgrid(x_, y_);

px0 = min(x_) + (max(x_) - min(x_)).*rand(1, 200);
py0 = min(y_) + (max(y_) - min(y_)).*rand(1, 200);

A0 = ones(1, 200);

vx = normrnd(5, 0.1, [1 200]); 
vy = normrnd(0, 0, [1 200]);
% vx = 3; vy = 0;
px(1, :) = px0; py(1, :) = py0;
A(1, :) = A0;

Nt = 64; dt = 1e-3;
lambda = 0.03;

% D(1, :) = sqrt(px(1, :).^2 + py(1, :).^2);
D(1, :) = zeros(1, 200);
z(1) = sum(A0 .* exp(1j .* 4*pi/lambda .* D(1, :))) + 0.001.^2 .* rand;


for i = 2:Nt
    px(i, :) = px(i - 1, :) + vx .* dt;
    py(i, :) = py(i - 1, :) + vy .* dt;
    
    inval_x = (px(i, :) > 30); inval_y = (py(i, :) > 30);
    [~, idx] = find(inval_x == 1); [~, idy] = find(inval_y == 1);
    Nx_inval = length(idx); Ny_inval = length(idy);
    
    if Nx_inval > 1
        px(i, idx) = min(x_) + (max(x_) - min(x_)).*rand(1, Nx_inval);
    end
    if Ny_inval > 1
        py(i, idy) = min(y_) + (max(y_) - min(y_)).*rand(1, Ny_inval);
    end
    
    D(i, :) = sqrt(px(i, :).^2 + py(i, :).^2);
    A(i, :) = A0;
    z(i) = sum(A(i, :) .* exp(-1j .* 4*pi/lambda .* D(i, :))) + 0.001.^2 .* rand;

end

ZFFT = 1./sqrt(Nt) .* fftshift(fft(z));

v_amb = lambda/(4 * dt);
dv = lambda/(4 * dt * (Nt - 1));
vel_axis = linspace(-v_amb, v_amb, Nt);

figure; plot(vel_axis, db(abs(ZFFT)));

PT = sum(abs(ZFFT).^2  .* dv);
mu = 1./PT .* sum(vel_axis .* abs(ZFFT).^2 .* dv);
sigma = sqrt(1./PT .* sum((vel_axis - mu).^2 .* abs(ZFFT).^2 .* dv));
