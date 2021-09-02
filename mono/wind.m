%% Wind field creation 
clear;
close all;

dr = 227;
[x, y] = meshgrid(-15e3:dr:15e3,  -15e3:dr:15e3);

%% velocity model ground truth

[u, v] = meshgrid(normrnd(400, 1, [size(x, 1) 1]), normrnd(0, 1, [size(x, 1) 1])); 

V = sqrt(abs(u).^2 + abs(v).^2);




% h = pcolor(x, y, V);
% set(h,'ZData',-1+zeros(size(V))); hold on;  shading flat; 
% quiver(x(1:e:end, 1:e:end), y(1:e:end, 1:e:end), u(1:e:end, 1:e:end), v(1:e:end, 1:e:end), 'w');

%% Reflectivity ground truth

Nx = size(x, 1); Ny = size(y, 2);

%% cloud1

x1 = 3000;
y1 = 1000;

x2 = -1000;
y2 = -900;

e = 0.8; 

a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
b = a*sqrt(1-e^2);
t = linspace(0,2*pi, 2000);
Xe = a*cos(t);
Ye = b*sin(t);
w = atan2(y2-y1,x2-x1);

xe = (x1+x2)/2 + Xe*cos(w) - Ye*sin(w);
ye = (y1+y2)/2 + Xe*sin(w) + Ye*cos(w);

for i = 1:length(xe)

    [~, idx(i)] = min(abs(x(1, :) - xe(i)));
    [~, idy(i)] = min(abs(y(:, 1) - ye(i)));
end

refx_ = zeros(size(x));
refx_(idy, idx) = 1;


figure;

h = pcolor(x, y, refx_);
set(h,'ZData',-1+zeros(size(refx_))); hold on;  shading flat; 
ed = 10;
quiver(x(1:ed:end, 1:ed:end), y(1:ed:end, 1:ed:end), u(1:ed:end, 1:ed:end), v(1:ed:end, 1:ed:end), 'w');


%% variation with time

Nt = 10;

for i = 1:Nt
    
    [~, idx1t] = min(abs(x(1, :) - x1));
    [~, idy1t] = min(abs(y(:, 1) - y1));
    
    x1 = x1 + u(1, idx1t) * i;
    y1 = y1 + v(idy1t, 1) * i;
    
    [~, idx2t] = min(abs(x(1, :) - x2));
    [~, idy2t] = min(abs(y(:, 1) - y2));
    
    x2 = x2 + u(1, idx2t) * i;
    y2 = y2 + v(idy2t, 1) * i;
    
    
    e = 0.8; 

    a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
    b = a*sqrt(1-e^2);
    t = linspace(0,2*pi, 2000);
    Xe = a*cos(t);
    Ye = b*sin(t);
    w = atan2(y2-y1,x2-x1);

    xe = (x1+x2)/2 + Xe*cos(w) - Ye*sin(w);
    ye = (y1+y2)/2 + Xe*sin(w) + Ye*cos(w);

    for k = 1:length(xe)

        [~, idx(k)] = min(abs(x(1, :) - xe(k)));
        [~, idy(k)] = min(abs(y(:, 1) - ye(k)));
    end

    refx_ = zeros(size(x));
    refx_(idy, idx) = 1;


    figure;

    h = pcolor(x, y, refx_);
    set(h,'ZData',-1+zeros(size(refx_))); hold on;  shading flat; 
    ed = 10;
    quiver(x(1:ed:end, 1:ed:end), y(1:ed:end, 1:ed:end), u(1:ed:end, 1:ed:end), v(1:ed:end, 1:ed:end), 'w');

end

%% radar

r = sqrt(x.^2 + y.^2);
phi = angle(x + 1j .* y);





