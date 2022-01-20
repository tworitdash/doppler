%% Vortex wind field

dph = 1.8 * pi/180;
R = 2e3;
Phi = 0*pi/180;
X0 = R .* cos(Phi);
Y0 = R .* sin(Phi);


r_ = eps:100:5e3;
phi_ = eps:dph:2*pi; 


[r, phi] = meshgrid(r_, phi_);


U = 5 .* ones(size(r));
V = 5 .* ones(size(r));

txt = [' Wind field '];

xl = 'x [km]';
yl = 'y [km]';
zl = 'U \hat x + V \hat y [m.s^{-1}]';
ex = 5;
ey = 5;
figure;    
quiverr(r.*cos(phi), r.*sin(phi), U, V, xl, yl, ex, ey, txt); 


x = r .* cos(phi); 
y = r .* sin(phi);
rx0 = 2000; ry0 = 2000;

xd = (x - X0)/rx0;
yd = (y - Y0)/ry0;


rd = sqrt(xd.^2 + yd.^2);
thetap = 5*pi/180;

theta = pi/2 - atan2(yd,rd) + thetap; sigma = 1;

A = 5 .* exp(-(rd - 1)/(2 .* sigma.^2));
u = A .* cos(theta); v = A .* sin(theta);

Ud = U + u; Vd = V + v;


txt = [' Wind field '];

xl = 'x [km]';
yl = 'y [km]';
zl = 'U \hat x + V \hat y [m.s^{-1}]';
ex = 10;
ey = 10;
figure;
quiverr(r.*cos(phi), r.*sin(phi), Ud, Vd, xl, yl, ex, ey, txt); 
