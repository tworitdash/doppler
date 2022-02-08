
%% Correlation in space and time
clear;
close all;

dd = 1;
dt = 15*60;

dx = 1e3; 
dy = 1e3;
Nx = 25;
Ny = 25;
% x_ = 1:dx:Nx*dx;
% y_ = 1:dy:Ny*dy;
x_ = linspace(0, 25e3, Nx); y_ = linspace(0, 25e3, Ny); 
Nxy = Nx * Ny;
% [x, y] = meshgrid(x_, y_);
x = repmat(x_, [1 Ny]);
y = repmat(y_, [1 Nx]);

X.x = x;
X.y = y;

% delta = sqrt((x - x(1)).^2 + (y - y(1)).^2);

% tau = eps:dt:5;
Nt_ = 16;

tau_v = 0:dt:Nt_*dt;

Nt = length(tau_v);


W.u = 1 * ones(1, Nxy);
W.v = 1 * ones(1, Nxy);
K = zeros(Nt, Nxy, Nxy);
Q = zeros(Nt, Nxy, Nxy);

for t = 1:Nt
    for i = 1:Nxy
        for k = 1:Nxy
            xi = X.x(i);
            yi = X.y(i);
            xk = X.x(k);
            yk = X.y(k);

            Xi = [xi yi];
            Xk = [xk yk];

            u = W.u(i);
            v = W.v(i);

            Wt = [u v] .* tau_v(t); 
%             Wt = [1 0];
            
            D = eye(2, 2);

            K(t, i, k) = 0.33 * exp(-(Xi - Xk - Wt) * inv(D) * (Xi - Xk - Wt).');

%             if i == k
%                 K(t, i, k) = 1;
%             else
%                 K(t, i, k) = 0;
%             end
            
            sigma = sqrt(1e-3);
            gamma = 3.33;
            Dik = sqrt(sum((Xi - Xk).^2));
            Q(t, i, k) = sigma^2 * exp(-Dik/gamma);
        end
    end
end

% figure; surface(1:Nxy, [1:Nxy].', squeeze(g(1, :, :))); shading flat;

%%
[R0] = RF0(x_, y_, 0.1, 0.1);
% R0 = rand(Nx, Ny);
R = zeros(Nt, Nx, Ny);
R(1, :, :) = R0./max(max(R0));

for t = 2:Nt
    R0_vec = reshape(squeeze(R(t - 1, :, :)), [Nxy 1]);
    qt =  mvnrnd(zeros(1, Nxy), squeeze(Q(t, :, :)));
    R1_vec = squeeze(K(t - 1, :, :)) * R0_vec;
    R(t, :, :) = reshape(R1_vec, Nx, Ny) + reshape(qt, Nx, Ny);
end

%% Plot fields
figure; h = surface(x_, y_, abs(squeeze(R(1, :, :)).')); shading interp; colorbar;
for t = 2:Nt
    h.CData = abs(squeeze(R(t, :, :))).';
    %figure; surface(x_, y_, abs(squeeze(R(t, :, :)).')); shading interp; colorbar;
    caption = sprintf('Frame #%d of %d, t = %.1f', t, Nt, tau_v(t));
	title(caption, 'FontSize', 15);
% 	drawnow;
	thisFrame = getframe(gca);
	% Write this frame out to a new video file.
%  	writeVideo(writerObj, thisFrame);
% 	myMovie(t) = thisFrame;
    
end

% movie(myMovie);
