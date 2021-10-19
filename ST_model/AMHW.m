%% Correlation in space and time
clear;
close all;

dd = 1;
dt = 1;

dx = 1; 
dy = 1;

x = 1:dx:75;
y = 1:dy:75;

% delta = sqrt((x - x(1)).^2 + (y - y(1)).^2);

% tau = eps:dt:5;

tau_v = eps:dt:32*dt;

% tau_v = 0;

bs = 5;
cs = 1;
bt = 3;
ct = 1;
theta = 0;

Theta = [bs cs bt ct theta];

vx = 4; vy = 0;


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
Nx = length(x); Ny = length(y);
Nxy = Nx * Ny;
Nt = length(tau_v);

DL_ = reshape(squeeze(DL(2, :, :, :, :)), [Nxy Nxy]);

DeltaLT = reshape(squeeze(DL), [Nt Nxy Nxy]);

figure; imagesc(1:Nxy, [1:Nxy].', DL_); shading flat;

for t = 1:length(tau_v)
    [K(t, :, :)] = AMHW_corr(Theta, squeeze(DeltaLT(t, :, :)), tau_v(t));
end

figure; surface(1:Nxy, [1:Nxy].', squeeze(K(1, :, :))); shading flat;



% for k = 1:Nxy
%        for m = 1:Nxy
%            deltaL = DeltaLT(:, k, m);
%            K = AMHW_corr(Theta, deltaL, tau_v.');
%            Md1 = varm(1, 4);
%            Est = estimate(Md1, [K]);
%            AR(k, m, :) = cell2mat(Est.AR);
%        end
% end

%% Initial random field 
% 
% [x_, y_]  = meshgrid(x, y);
% 
% R0 = rand(size(x_));
% 
% R0 = fft2(R0);
% 
% 
% fx = 0.05;
% fy = 0.02; 
% 
% nx1 = 1 + floor(0.9999.*fx*length(x));
% nx2 = length(x) - nx1 + 1; 
% 
% ny1 = 1 + floor(0.9999.*fy*length(y));
% ny2 = length(y) - ny1 + 1;
% 
% 
% R0(nx1:nx2, :) = 0;
% R0(:, ny1:ny2) = 0;
% 
% R0 = ifft2(R0);
% 
% Rt = zeros(Nt, length(x), length(y));
% R(1, :, :) = R0;
% 
% for t = 1:Nt
%     for i = 1:length(x)
%         for l = 1:length(y)
%             AR_ = AR()
%             if (t > 3)
%                 R(t, :, :) = squeeze(sum(R(t-4:t, :, :) .* AR, 1));
%             elseif (t > 1) && (t < 3)
%                 R(t, :, :) = squeeze(sum(R(1:t, :, :) .* AR(:,:, 1:t)));
%             else
%                 R(t, :, :) = squeeze(sum(R(1:t, :, :) .* AR(:,:, 1:t)));
%             end
%         end
%     end
% end 


%% Cholskey on time axis
A = zeros(Nt, Nxy, Nxy);
RF = zeros(Nt, Nxy);
RF_rv = zeros(Nt, Nx, Ny);

[R1] = RF0(x, y, 0.07, 0.07);

RF(1, :) = reshape(R1, [1 Nxy]);
% A(1, :, :) = chol((squeeze(K(1, :, :))));
RF_rv(1, :, :) = R1;

figure; surface(x, y, squeeze(abs(RF_rv(1 , :, :)))); shading flat; 

for t = 2:Nt
%     U = nearestSPD(squeeze(K(t, :, :)));
%     A(t, :, :) = chol((U));
    RF(t, :) = ((squeeze(K(t - 1, :, :)))) * RF(t-1, :).';
    RF_rv(t, :, :) = reshape(RF(t, :), [Nx Ny]);
end

%%
for t = 1:Nt
%     figure; surface(1:Nxy, 1:Nxy, squeeze(K(t, :, :))); shading flat; 
    figure; surface(x, y, squeeze(abs(RF_rv(t, :, :)))); shading flat; 
end