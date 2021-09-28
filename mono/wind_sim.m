%% Wind field creation 
clear;
close all;

dr = 227;
[x, y] = meshgrid(-15e3:dr:15e3,  -15e3:dr:15e3);

%% velocity model ground truth

[u, v] = meshgrid(normrnd(40, 1, [size(x, 1) 1]), normrnd(40, 10, [size(x, 1) 1])); 

V = sqrt(abs(u).^2 + abs(v).^2);




% h = pcolor(x, y, V);
% set(h,'ZData',-1+zeros(size(V))); hold on;  shading flat; 
% quiver(x(1:e:end, 1:e:end), y(1:e:end, 1:e:end), u(1:e:end, 1:e:end), v(1:e:end, 1:e:end), 'w');

%% Reflectivity ground truth

Nx = size(x, 1); Ny = size(y, 2);

%% cloud1

refx_ = zeros(size(x));
% input the vertices of all ellipses 

x1_ = [1e3 -3e3 -2e3 7e3 -2e3 ];
y1_ = [10e3 10e3 -8e3 -12e3 0 ];

x2_ = [4e3 -10e3 -8e3 12e3 -5e3 ];
y2_ = [10e3 10e3 -8e3 -12e3 0 ];

for m = 1:length(x1_)

e = 0.8; 

a = 1/2*sqrt((x2_(m)-x1_(m))^2+(y2_(m)-y1_(m))^2);
b = a*sqrt(1-e^2);
t = linspace(0,2*pi, 2000);
Xe = a*cos(t);
Ye = b*sin(t);
w = atan2(y2_(m)-y1_(m),x2_(m)-x1_(m));

xe = (x1_(m)+x2_(m))/2 + Xe*cos(w) - Ye*sin(w);
ye = (y1_(m)+y2_(m))/2 + Xe*sin(w) + Ye*cos(w);

for i = 1:length(xe)

    [~, idx(i)] = min(abs(x(1, :) - xe(i)));
    [~, idy(i)] = min(abs(y(:, 1) - ye(i)));
end


refx_(idy, idx) = 1;

end

figure;

h = pcolor(x, y, refx_);
set(h,'ZData',-1+zeros(size(refx_))); hold on;  shading flat; 
ed = 10;
quiver(x(1:ed:end, 1:ed:end), y(1:ed:end, 1:ed:end), u(1:ed:end, 1:ed:end), v(1:ed:end, 1:ed:end), 'w');


%% variation with time

Nt = 10;

for i = 1:Nt
    for m = 1:length(x1_)
        [~, idx1t] = min(abs(x(1, :) - x1_(m)));
        [~, idy1t] = min(abs(y(:, 1) - y1_(m)));

        x1_(m) = x1_(m) + u(1, idx1t) * i;
        y1_(m) = y1_(m) + v(idy1t, 1) * i;

        [~, idx2t] = min(abs(x(1, :) - x2_(m)));
        [~, idy2t] = min(abs(y(:, 1) - y2_(m)));

        x2_(m) = x2_(m) + u(1, idx2t) * i;
        y2_(m) = y2_(m) + v(idy2t, 1) * i;


        e = 0.8; 

        a = 1/2*sqrt((x2_(m)-x1_(m))^2+(y2_(m)-y1_(m))^2);
        b = a*sqrt(1-e^2);
        t = linspace(0,2*pi, 2000);
        Xe = a*cos(t);
        Ye = b*sin(t);
        w = atan2(y2_(m)-y1_(m),x2_(m)-x1_(m));
   
    
        xe = (x1_(m)+x2_(m))/2 + Xe*cos(w) - Ye*sin(w);
        ye = (y1_(m)+y2_(m))/2 + Xe*sin(w) + Ye*cos(w);

        for k = 1:length(xe)

            [~, idx(m, k)] = min(abs(x(1, :) - xe(k)));
            [~, idy(m, k)] = min(abs(y(:, 1) - ye(k)));
        end
        
        
       
    end
    
    
    idx_ = reshape(idx, [1, 5 * 2000]);
    idy_ = reshape(idy, [1, 5 * 2000]);
    
    
    refx_ = zeros(size(x));
    refx_(idy_, idx_) = 1;
    
    figure;

    h = pcolor(x, y, refx_);
    set(h,'ZData',-1+zeros(size(refx_))); hold on;  shading flat; 
    ed = 10;
    quiver(x(1:ed:end, 1:ed:end), y(1:ed:end, 1:ed:end), u(1:ed:end, 1:ed:end), v(1:ed:end, 1:ed:end), 'w');

end

%% radar

r = sqrt(x.^2 + y.^2);
phi = angle(x + 1j .* y);






