clear;
close all;

%% Random Fields in 2D


dx = 227;
dy = 227;
x_ = -15e3:dx:15e3;
y_ = -15e3:dy:15e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

% Min = 0.0001;
Max = 1;

R0 = Max .* (rand(size(x)) + 1j .* rand(size(x)));

% R0 = (rand(size(x))./sqrt(2) + 1j .* rand(size(x))./sqrt(2));


figure; imagesc(x_, y_, (abs(R0))); colormap('jet')

shading flat; colorbar; axis equal tight;

%% GENERATE RANDOM FUNCTION 
SNR_db = 40;
SNR = 10^(SNR_db/10);

f = 0.02;

% rand('seed', 0);

nx1 = 1 + floor(0.9999.*f*length(x_));
nx2 = length(x_) - nx1 + 1; 

ny1 = 1 + floor(0.9999.*f*length(y_));
ny2 = length(y_) - ny1 + 1; 

F = 10^4 .* (rand(length(x), length(y)) + 1j .* rand(length(x), length(y)) - 0.5);
F = fft2(F);
figure; surface(x_.*1e-3, y_.*1e-3, db(abs(F).')./2); colorbar; shading interp;

F(nx1:nx2, :) = 0;
F(:, ny1:ny2) = 0;

F             = (ifft2(F));
N             = sqrt(max(abs(F)).^2./SNR) .* (rand(length(x), length(y)) + 1j .* rand(length(x), length(y)));
F             = F + N;

pcolor(x_.*1e-3, y_.*1e-3, db(abs(F).')./2); colorbar; 

% pcolor(x_.*1e-3, y_.*1e-3, db(abs(N).')./2); colorbar; 
shading interp; axis equal tight; 

xlim([-15 15]);
ylim([-15 15]);

T = [-15:5:15];
L = {};
for m = 1:length(T)
    if T(m) == 0
        L{m} = '0';
    else
        L{m} = num2str(T(m), '%3.1f');
    end
end

set(gca, 'XTick', T, 'XTickLabel', L);
set(gca, 'YTick', T, 'YTickLabel', L);

set(gca, 'FontSize', 18, 'LineWidth', 4);

title('Relectivity $\eta(x, y)$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);

%% VELOCITY FIELDS


