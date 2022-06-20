clear;
close all;

%% Random Fields in 2D


dx = 25;
dy = 25;
x_ = -dx*100:dx:dx*100;
% y_ = -5e3:dy:5e3;
y_ = -dy*100:dy:dy*100;
% x_ = -100:dx:100;
% y_ = -100:dy:100;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle
lambdax = 10;
lambday = 10;

fx = 0.1;
fy = 0.1;

[E0] = Eta_0(x_, y_, x, y, dx, dy, lambdax, lambday, fx, fy);

E0_norm = E0./max(max(E0));

Energy0 = sum(sum(abs(E0_norm).^2 .* dx .* dy));

txt = ['Initial reflectivity field'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field Intensity';


surplot(x*1e-3, y*1e-3, abs(E0_norm), xl, yl, zl, txt);
% quiver(x*1e-3, y*1e-3, E0_norm, 0.2 * E0_norm);


%% R1 is the cloud reflectivity field


U = normrnd(5, 0.1, size(x));
V = normrnd(0, 0, size(y));


dt = 1; % Volume Scan Time
Nt = 10; % Number of time steps

[Rt, Eta] = RFTV3(x_, y_, x, y, dx, dy, E0, U, V, dt, Nt, Energy0);

t = eps:dt:Nt*dt;


close all;

Enable_movie = input('Do you want yo plot the movie of ground trurh? [0]: , ');

idx = find((x_<1.1e2)&(x_>-1.1e2)); idy = find((y_ < 1.1e2)&(y_>-1.1e2)); 
xd_ = x_(idx); yd_ = y_(idy);

[xd, yd] = meshgrid(xd_, yd_);

if Enable_movie == 1
figure;
h = surface(xd.*1e-3, yd.*1e-3, abs(squeeze(Rt(1, idy, idx))));


    for k = 2:Nt

        
        h.CData = abs(squeeze(Rt(k, idy, idx)));
        
        shading flat; colormap('jet'); colorbar; % caxis([0.9 1]);
        pause(0.1);
        caption = sprintf('Frame #%d of %d, t = %.1f', k, Nt, t(k));
        title(caption, 'FontSize', 15); 

        
    end
end
%% RADAR Filter

% Rt_req = Rt(:, :, :);
% R = 1e2;
% dr = 20; dph = 1.8*pi/180;
% Omega_rpm = 1;
% Omega = Omega_rpm * 2 * pi / 60;
% lambda = 0.03;
% phi = eps:dph:3*dph;
% r = eps:dr:R;
% 
% [Rr, N_sweep] = Radar_Filter_V2(r, phi, xd_, yd_, xd, yd, dx, dy, Rt_req, Nt, dt, Omega, dph, lambda);
% 
% [PT, Mu, Sigma] = Ret(Rr, r, phi, N_sweep, dt, lambda); 
% 
% %% Plot radar image with time
% 
% close all;
% 
% [rr, phir] = meshgrid(r, phi);
% 
% xr = rr .* cos(phir);
% yr = rr .* sin(phir);
% 
% txt = ['Total Power'];
% % 
% xl = 'x [km]';
% yl = 'y [km]';
% zl = '\mu [m.sec^{-1}]';
% 
% 
% surplot(xr*1e-3, yr*1e-3, db(PT.'), xl, yl, zl, txt);
% 
% txt = ['Radial Mean Doppler Velocity [m.sec^{-1}]'];
% % 
% xl = 'x [km]';
% yl = 'y [km]';
% zl = '\mu [m.sec^{-1}]';
% 
% 
% surplot(xr*1e-3, yr*1e-3, Mu.', xl, yl, zl, txt);
% 
% txt = ['Radial Doppler Spectrum Width [m.sec^{-1}]'];
% % 
% xl = 'x [km]';
% yl = 'y [km]';
% zl = '\sigma [m.sec^{-1}]';
% 
% 
% surplot(xr*1e-3, yr*1e-3, Sigma.', xl, yl, zl, txt);