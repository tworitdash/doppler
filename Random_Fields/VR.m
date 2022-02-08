clear;
close all;

% dx = 1000; 
% dy = 1000;

% x_ = -25e3:dx:25e3;
% y_ = -25e3:dy:25e3;

dx = 1e-3; 
dy = 1e-3;

BW_deg = 1.8;
BW = BW_deg .* pi/180;
X_end = 5;

x_ = eps:dx:X_end;
y_ = eps:dy:X_end*tan(BW);

[x, y] = meshgrid(x_, y_);

% [E0] = Eta_0(x_, y_, x, y, dx, dy, 2e3, 2e3, 0.07, 0.07);

[E0] = Eta_0(x_, y_, x, y, dx, dy, 0.1, 0.1, 0.1, 0.1);
E0_norm = E0./max(max(E0));

Energy0 = sum(sum(abs(E0_norm).^2 .* dx .* dy));

txt = ['Initial reflectivity field'];
% 
% xl = 'x [km]';
% yl = 'y [km]';
% zl = 'Random Field Intensity';
% surplot(x*1e-3, y*1e-3, abs(E0_norm), xl, yl, zl, txt);

xl = 'x [m]';
yl = 'y [m]';
zl = 'Random Field Intensity';

surplot(x, y, abs(E0_norm), xl, yl, zl, txt);
Nx = length(x_);
Ny = length(y_);

Nxy = Nx * Ny;
 
D = eye(2, 2);
alpha = 0.33;
% dt = 15*60;
dt = 1e-3;
W = [0 1]*dt;
sigma_s = 1e-1;
Gamma = 0.33;
theta = 0.5;
for i = 1:Nxy
    for k = 1:Nxy
                x_1 = x(i); x_2 = x(k);
                y_1 = y(i); y_2 = y(k);

                X_1 = [x_1 y_1]; X_2 = [x_2 y_2];
%                 Corr = exp(-(X_1 - W - X_2) * inv(D) * (X_1 - W - X_2).'./2e3);
                Corr = exp(-(X_1 - W - X_2) * inv(D) * (X_1 - W - X_2).');
                H(i, k) = alpha .* Corr;
                Dik = ((X_1 - W - X_2) * (X_1 - W - X_2).');
                Q(i, k) = sigma_s.^2 .* exp(-(Dik)./(2 * Gamma.^2));
    end
end

Rt(1, :, :) = E0;
Nt = 8;
t = eps:dt:Nt*dt;


for ti = 1:Nt - 1
   utb =  reshape(Rt(ti, :, :), [Nxy 1]);
   uta = H * utb + sigma_s.^2 * mvnrnd(zeros(Nxy, 1), Q).';
   Rt(ti+1, :, :) = reshape(uta, [Nx Ny]);
end
%%

figure;
h = surface(x.*1e-3, y.*1e-3, abs(squeeze(Rt(1, :, :))).');


    for m = 2:Nt

        
        h.CData = abs(squeeze(Rt(m, :, :))).';
        
%         shading flat; 
        colormap('jet'); colorbar; % caxis([0.9 1]);
        pause(0.1);
        caption = sprintf('Frame #%d of %d, t = %.1f', m, Nt, t(m));
        title(caption, 'FontSize', 15); 

        
    end

%% 
mp = [1 2 4 6 8];

for i = 1:length(mp)
    
txt = ['Reflectivity at time t = ', num2str(t(mp(i))), ' [sec]'];
% 
xl = 'x [km]';
yl = 'y [km]';
zl = 'Random Field';

surplot(x*1e-3, y*1e-3, abs(squeeze(Rt(mp(i), :, :))).', xl, yl, zl, txt);

end
