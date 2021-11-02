clear;
close all;

%% Random Fields in 2D


dx = 100;
dy = 100;
x_ = -25e3:dx:25e3;
y_ = -25e3:dy:25e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

[E0] = Eta_0(x_, y_, x, y, dx, dy, 1e3, 0.1);

E0_norm = E0./max(max(E0));

txt = ['Initial random field'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field Intensity';


surplot(x, y, abs(E0_norm), xl, yl, zl, txt);


%% R1 is the cloud reflectivity field
U = normrnd(20, 0.1, size(x));
V = normrnd(20, 0.1, size(x));


dt = 1; % Volume Scan Time
Nt = 100; % Number of time steps

[Rt] = RFTV3(x_, y_, x, y, dx, dy, E0, U, V, dt, Nt);

t = eps:dt:Nt*dt;


close all;

Enable_movie = input('Do you want yo plot the movie of ground trurh? [0]: , ');

idx = find((x_<15.1e3)&(x_>-15.1e3)); idy = find((y_ < 15.1e3)&(y_>-15.1e3)); 
xd_ = x_(idx); yd_ = y_(idy);

[xd, yd] = meshgrid(xd_, yd_);

if Enable_movie == 1
figure;
h = surface(xd, yd, abs(squeeze(Rt(1, idy, idx))).');


    for k = 2:Nt

        
        h.CData = abs(squeeze(Rt(k, idy, idx))).';
        
        shading flat; colormap('jet'); colorbar; % caxis([0.9 1]);
        pause(0.1);
        caption = sprintf('Frame #%d of %d, t = %.1f', k, Nt, t(k));
        title(caption, 'FontSize', 15); 

        
    end
end
%% RADAR Filter

lambda = 0.03;
R = 5e3;
dr = 100; dph = 2*pi/180; 

phi = eps:dph:pi/2;
r = eps:dr:R;

[Rr] = Radar_Filter(r, phi, xd_, yd_, xd, yd, dx, dy, Rt(:, idy, idx), Nt);

%% Plot radar image with time 
close all;

[rr, phir] = meshgrid(r, phi);

xr = rr .* cos(phir);
yr = rr .* sin(phir);

Enable_movie = input('Do you want yo plot the movie of the radar image? [0]: , ');

if Enable_movie == 1
figure;
h = surface(xr, yr, abs(squeeze(Rr(1, :, :))).');


    for k = 2:Nt

        
        h.CData = abs(squeeze(Rr(k, :, :))).';
        
        shading flat; colormap('jet'); colorbar; % caxis([0.9 1]);
        pause(0.1);
        caption = sprintf('Frame #%d of %d, t = %.1f', k, Nt, t(k));
        title(caption, 'FontSize', 15); 

        
    end
end



% figure; surface(xr, yr, (abs(squeeze(Rr(1, :, :))).')./2); colorbar; colormap('jet'); shading flat;
%% Random Field and Radar Filter comparison

txt = ['Random field'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field';

m = 5;
surplot(x*1e-3, y*1e-3, abs(squeeze(Rt(m, :, :))).', xl, yl, zl, txt);


txt = ['Random field after Radar Filter'];
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field';

surplot(xr*1e-3, yr*1e-3, abs(squeeze(Rr(m, :, :))).', xl, yl, zl, txt);
