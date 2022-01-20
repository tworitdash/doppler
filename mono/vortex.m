x = 0:0.1:10;
z = 0:0.1:5;

[x, z] = meshgrid(x, z);

t = 0;

x0 = 4; 
z0 = 0.5;

rx0 = 2; 
rz0 = 2;

xst = (x - x0)/rx0;
zst = (z - z0)/rz0;

r = sqrt(xst.^2 + zst.^2);
theta_p = pi/2;

theta = pi/2 - atan2(zst, r) + theta_p;
sigma = 1;
A_0 = 5;
A = A_0 .* exp(-(r - 1)./(2 .* sigma^2));

U = 0 .* ones(size(x)); V = 0 .* ones(size(x));

u = U +  A .* cos(theta); v = V - A .* sin(theta);
e = 1;
figure; quiver(x(1:e:end, 1:e:end), z(1:e:end, 1:e:end), u(1:e:end, 1:e:end), v(1:e:end, 1:e:end)); 