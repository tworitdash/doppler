clear;

%% radar specification

dr = 50;
dph_deg = 1.8;
dph = dph_deg*pi/180;

r = eps:dr:5e3;
phi = eps:dph:2*pi;
theta = 0;

%% Analytical Reflectivity 

Lambda = 1e3;

[r_, phi_] = meshgrid(r, phi);
dT = 10;
n_rot = 100;
t = eps:dT:dT*n_rot;

phi_0 = 0;

U = 1; V = 1;

Vr = U .* cos(theta) .* cos(phi_) + V .* cos(theta) .* sin(phi_);

Vph = U .* cos(theta) .* sin(phi_) - V .* cos(theta) .* cos(phi_);

figure;

for i = 1:length(t)

Eta(i, :, :) = sin(2*pi/Lambda * (r_ - Vr .* t(i))) .* cos(phi_ - Vph./Vr .* log(r_));

surface(r_ .* cos(theta) .* cos(phi_), r_ .* cos(theta) .* sin(phi_), abs(squeeze(Eta(i, :, :)))); shading flat; pause(1);

end

