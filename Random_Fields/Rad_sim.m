%% Ground truth
clear;
dx = 10;
dy = 10;

x_ = eps:dx:500;
y_ = eps:dy:500;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

[E0] = Eta_0(x_, y_, x, y, dx, dy, 1e2, 1e2, 0.1, 0.1);

E0_norm = E0./max(max(E0));

Energy0 = sum(sum(abs(E0_norm).^2 .* dx .* dy));

txt = ['Initial reflectivity field'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field Intensity';


surplot(x*1e-3, y*1e-3, abs(E0_norm), xl, yl, zl, txt);

%% 
Nx = length(x_);
Ny = length(y_);

Nxy = Nx * Ny;
 
D = eye(2, 2);
alpha = 0.33;
dt = 1e-3;
W = [9 0]*dt;
sigma_s = 1e-1;
Gamma = 0.33;
theta = 0.5;


for i = 1:Nxy
    for k = 1:Nxy
                x_1 = x(i); x_2 = x(k);
                y_1 = y(i); y_2 = y(k);

                X_1 = [x_1 y_1]; X_2 = [x_2 y_2];
                Corr = alpha .* exp(-(X_1 - W - X_2) * inv(D) * (X_1 - W - X_2).');
                H(i, k) = alpha .* Corr;
                Dik = ((X_1 - W - X_2) * (X_1 - W - X_2).');
                Q(i, k) = sigma_s.^2 .* exp(-(Dik)./(2 * Gamma.^2));
    end
end

Rt(1, :, :) = E0;
Nt = 2048;
t = eps:dt:Nt*dt;


for ti = 1:Nt - 1
   utb =  reshape(Rt(ti, :, :), [Nxy 1]);
   uta = H * utb; % + sigma_s.^2 * mvnrnd(zeros(Nxy, 1), Q).';
   Rt(ti+1, :, :) = reshape(uta, [Ny Nx])./max(abs(uta));
end

figure;
h = surface(x.*1e-3, y.*1e-3, abs(squeeze(Rt(1, :, :))));


    for m = 2:Nt

        
        h.CData = abs(squeeze(Rt(m, :, :)));
        
%         shading flat; 
        colormap('jet'); colorbar; % caxis([0.9 1]);
        pause(0.1);
        caption = sprintf('Frame #%d of %d, t = %.1f', m, Nt, t(m));
        title(caption, 'FontSize', 15); 

        
    end
    
 %% Radar filter
 
 dr = 10; 
 dph = 1.8*pi/180;
 ph_ = eps:dph:2*dph;
 r_ = eps:dr:500;
 
 [r, ph] = meshgrid(r_, ph_);
 
 figure; surface(r.*cos(ph), r.*sin(ph), ones(length(r_), length(ph_)).'); shading flat;
 