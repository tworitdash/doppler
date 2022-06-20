clear;
close all;

%% Random Fields in 2D


dx = 10;
dy = 10;
x_ = -5e3:dx:5e3;
y_ = -5e3:dy:5e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

[E0] = Eta_0(x_, y_, x, y, dx, dy, 1e3, 0.5, 0.5, 0.5);

E0_norm = E0./max(max(E0));

Energy0 = sum(sum(abs(E0_norm).^2 .* dx .* dy));

txt = ['Initial reflectivity field'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field Intensity';


surplot(x*1e-3, y*1e-3, abs(E0_norm).', xl, yl, zl, txt);


%% R1 is the cloud reflectivity field
U = normrnd(20, 0.1, size(x));
V = normrnd(20, 0.1, size(x));


dt = 1e-3; % Volume Scan Time
Nt = 128; % Number of time steps

[Rt, Eta] = RFTV3(x_, y_, x, y, dx, dy, E0, U, V, dt, Nt, Energy0);

t = eps:dt:Nt*dt;


close all;

Enable_movie = input('Do you want yo plot the movie of ground trurh? [0]: , ');

idx = find((x_<1.1e3)&(x_>-1.1e3)); idy = find((y_ < 1.1e3)&(y_>-1.1e3)); 
xd_ = x_(idx); yd_ = y_(idy);

[xd, yd] = meshgrid(xd_, yd_);

if Enable_movie == 1
figure;
h = surface(xd.*1e-3, yd.*1e-3, abs(squeeze(Rt(1, idy, idx))).');


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
PRT = 1e-3;
R = 1e3;
dr = 100; dph = 2*pi/180; 

phi = eps:dph:pi/2;
r = eps:dr:R;

Omega_rpm = 60;
Omega = Omega_rpm/(2*pi/60);

[Rr] = Radar_Filter_V2(r, phi, xd_, yd_, xd, yd, dx, dy, Rt(:, idy, idx), Nt, PRT, Omega, dph, lambda);
                     

%% Plot radar image with time 
close all;

[rr, phir] = meshgrid(r, phi);

xr = rr .* cos(phir);
yr = rr .* sin(phir);

Enable_movie = input('Do you want yo plot the movie of the radar image? [0]: , ');

if Enable_movie == 1
figure;
h = surface(xr.*1e-3, yr.*1e-3, abs(squeeze(Rr(1, :, :))).');


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
mp = [2 25 50 75 99];

for i = 1:length(mp)
    
txt = ['Reflectivity at time t = ', num2str(t(mp(i))), ' [sec]'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field';

surplot(xd*1e-3, yd*1e-3, abs(squeeze(Rt(mp(i), idy, idx))).', xl, yl, zl, txt);


%% 
txt = ['Random field after Radar Filter'];
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field';

surplot(xr*1e-3, yr*1e-3, abs(squeeze(Rr(mp(i), :, :))).', xl, yl, zl, txt);
end
