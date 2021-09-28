%% 
clear; 
close all;

lambda = 0.03;
%% Radar forward model 

th = linspace(2.5, 60, 15) * pi/180;
% ph = linspace(eps, 360, 360) * pi/180;
ph = linspace(eps, 1.88, 1024) * pi/180;

dR = 227;
R_max = 15e3;

r = eps:dR:R_max;

[R, phi, theta] = ndgrid(r, ph, th);

%% 

M0_th = Gauss_PDF(1, 10*pi/180, 2*pi/180, th) ...
    + Gauss_PDF(1, 25*pi/180, 2*pi/180, th);

% 
figure; 
plot(th*180/pi, db(M0_th)/2); grid on;


%% Sectors in azimuth with reflectivity

M0_ph_ini = zeros(size(ph));

% M01_ph = Gauss_PDF(10^(10/10), 10*pi/180, 2*pi/180, ph1) + ...
%     Gauss_PDF(10^(30/10), 25*pi/180, 2*pi/180, ph1);
% 
% M02_ph = Gauss_PDF(10^(10/10), 100*pi/180, 2*pi/180, ph2) + ...
%     Gauss_PDF(10^(30/10), 110*pi/180, 2*pi/180, ph2);
% 
% M03_ph = Gauss_PDF(10^(10/10), 235*pi/180, 2*pi/180, ph3) + ...
%     Gauss_PDF(10^(30/10), 250*pi/180, 2*pi/180, ph3);

M01_ph = (Antenna_BEAM(20 .* lambda, lambda, ph, 0.9*pi/180)).^2;
M01_ph = M01_ph./max(M01_ph);


M0_ph_ini = M01_ph; 

figure; plot(ph*180/pi, db(M01_ph));

% M0_ph_ini(1:30) = M01_ph;
% M0_ph_ini(91:120) = M02_ph;
% M0_ph_ini(225:260) = M03_ph;

%% 

M0_th_ph = M0_ph_ini.' * M0_th;


%% Range information 
% 
% [~, idxmin] = min(abs(r - 10e3));
% 
% [~, idxmax] = min(abs(r - 14e3));
% 
% [~, idxmin1] = min(abs(r - 1e3));
% 
% [~, idxmax1] = min(abs(r - 3e3));
% 
% r1 = r(idxmin:idxmax);
% r2 = r(idxmin1:idxmax1);

sigma_r = 0.35 .* dR;

M01_r = Gauss_PDF(1, 12e3, sigma_r, r) + Gauss_PDF(1, 2e3, sigma_r, r);

for i = 1:length(M01_r)
    M0(i, :, :) = M01_r(i) .* M0_th_ph;
end
%%

figure; surface(r, ph*180/pi, db(squeeze(M0(:, :, 1)).')); shading flat; colormap('jet'); colorbar;

%% M1 dependence

u_th = linspace(-5, 5, length(th));
v_th = linspace(0, 0, length(th));
w_th = linspace(0, 0, length(th));

% 
% for l = 1:length(M0_ph_ini)
%     if db(M0_ph_ini(l)) > -20
%         M0_ph_mask(l) = M0_ph_ini(l);
%     else
%         M0_ph_mask(l) = 0;
%     end
% end
% 
% for l = 1:length(M01_r)
%     if ( (l > idxmin) && (l < idxmax) ) || ( (l > idxmin1) && (l < idxmax1) )
%         M0_r_mask(l) = 1;
%     else
%         M0_r_mask(l) = 0;
%     end
% end

for k = 1:length(th)
    u_th_ph(:, k) = normrnd(u_th(k), 1, [1 length(ph)]) .* M0_ph_ini./max(M0_ph_ini);
    v_th_ph(:, k) = normrnd(v_th(k), 0, [1 length(ph)]) .* M0_ph_ini./max(M0_ph_ini);
    w_th_ph(:, k) = normrnd(w_th(k), 0.1, [1 length(ph)]) .* M0_ph_ini./max(M0_ph_ini);
end


for i = 1:length(M01_r)
    u(i, :, :) = M01_r(i)./max(M01_r) .* u_th_ph;
    v(i, :, :) = M01_r(i)./max(M01_r) .* u_th_ph;
    w(i, :, :) = M01_r(i)./max(M01_r) .* w_th_ph;
end

vr = u .* cos(theta) .* cos(phi) + v .* cos(theta) .* sin(phi) + w .* sin(theta);

% 
% figure; plot(ph * 180/pi, squeeze(vr(round((idxmin+idxmax)/2), :, 5))); grid on;


figure; quiver3(R .* cos(theta) .* cos(phi), R .* cos(theta) .* sin(phi), R .* sin(theta), u, v, w)

%% Constant elevation plots
ei = 3;

x_ei = squeeze(R(:, :, ei)) .* squeeze(cos(theta(:, :, ei))) .* squeeze(cos(phi(:, :, ei)));
y_ei = squeeze(R(:, :, ei)) .* squeeze(cos(theta(:, :, ei))) .* squeeze(sin(phi(:, :, ei)));
u_ei = squeeze(u(:, :, ei));
v_ei = squeeze(v(:, :, ei));
vr_ei = squeeze(vr(:, :, ei));
M0_ei = squeeze(M0(:, :, ei));

R_ei = squeeze(R(:, :, ei));
phi_ei = squeeze(phi(:, :, ei));

figure; quiver(x_ei, y_ei, u_ei, v_ei);

%% Surface Plots
% 
% figure; polarPcolor(r, ph * 180/pi, vr_ei); shading flat; colormap('jet'); colorbar;
% figure; polarPcolor(r, ph * 180/pi, db(M0_ei)); shading flat; colormap('jet'); colorbar;


figure; surface(x_ei, y_ei, vr_ei); shading flat; colormap('jet'); colorbar; grid on; view(2);
figure; surface(x_ei, y_ei, db(M0_ei)); shading flat; colormap('jet'); colorbar; grid on; view(2);
 
%% 