
%% Correlation in space and time
clear;
close all;

dd = 1;
dt = 1;

dx = 1; 
dy = 1;
Nx = 32;
Ny = 32;
x = 1:dx:Nx;
y = 1:dy:Ny;

% delta = sqrt((x - x(1)).^2 + (y - y(1)).^2);

% tau = eps:dt:5;
Nt_ = 16;

tau_v = 0:dt:Nt_*dt;

Nt = length(tau_v);

% tau_v = 0;

bs = 10;
cs = 1;
bt = 30;
ct = 1;
theta = 1;

Theta = [bs cs bt ct theta];

vx = 0.5; vy = 0.5;


for t = 1:length(tau_v)
    for k = 1:length(x)
       for m = 1:length(y)
            for l = 1:length(x)
                for n = 1:length(y)
                    DL(t, k, m, l, n) = sqrt((x(k) - tau_v(t) .* vx - x(l)).^2 ...
                        + (y(m) - tau_v(t) .* vy - y(n)).^2);
                end
            end
        end
    end
end




Nxy = Nx * Ny;

DL_ = reshape(squeeze(DL(2, :, :, :, :)), [Nxy Nxy]);

DeltaLT = reshape(squeeze(DL), [Nt Nxy Nxy]);

% figure; imagesc(1:Nxy, [1:Nxy].', DL_); shading flat;
% figure; imagesc(1:Nxy, [1:Nxy].', squeeze(DeltaLT(2, :, :))); shading flat;


alpha = 0.33;


for t = 1:length(tau_v)
%     [K(t, :, :)] = AMHW_corr(Theta, squeeze(DeltaLT(t, :, :)), tau_v(t));
   
    K(t, :, :) = alpha * exp(- squeeze(DeltaLT(t, :, :)).^2);
    gamma = 3.33;
    sigma = 1e-3;
    Q(t, :, :) = sigma^2 * exp(-DeltaLT(t, :, :)/gamma);
end
%%
% figure; l = surface(1:Nxy, [1:Nxy].', abs(squeeze(K(1, :, :)).')); shading interp; colorbar; colormap('jet');
% 
% for t = 2:Nt
%     l.CData = abs(squeeze(K(t, :, :)).'); 
%     caption = sprintf('Frame #%d of %d, t = %.1f', t, Nt, tau_v(t));
% 	title(caption, 'FontSize', 15);
% % 	drawnow;
% 	thisFrame = getframe(gca);
% 	% Write this frame out to a new video file.
% %  	writeVideo(writerObj, thisFrame);
% % 	myMovie(t) = thisFrame;
%     
% end

%% Initial field

[R0] = RF0(x, y, 0.1, 0.1);
% R0 = rand(Nx, Ny);
R = zeros(Nt, Nx, Ny);
R(1, :, :) = R0./max(max(R0));

for t = 2:Nt
%     qt =  mvnrnd(zeros(1, Nxy), squeeze(Q(t, :, :)));
    R0_vec = reshape(squeeze(R(t - 1, :, :)), [Nxy 1]);
    R1_vec = squeeze(K(t - 1, :, :)) * R0_vec;
    R1_vec = R1_vec./max(max(R1_vec));
    R(t, :, :) = reshape(R1_vec, Nx, Ny); % + reshape(qt, Nx, Ny);
end

%% Plot fields
figure; h = surface(x, y, abs(squeeze(R(1, :, :)).')); shading interp; colorbar; colormap('jet');
for t = 2:Nt
    h.CData = abs(squeeze(R(t, :, :)).'); 
%     contourf(x, y, abs(squeeze(R(t, :, :)).')); shading interp; colorbar; colormap('jet');
    caption = sprintf('Frame #%d of %d, t = %.1f', t, Nt, tau_v(t)); caxis([0.5 1]);
	title(caption, 'FontSize', 15);
% 	drawnow;
	thisFrame = getframe(gca);
	% Write this frame out to a new video file.
%  	writeVideo(writerObj, thisFrame);
% 	myMovie(t) = thisFrame;
    
end

% movie(myMovie);