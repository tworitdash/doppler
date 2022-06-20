clear;
close all;

%% Random Fields in 2D


dx = 100;
dy = 100;
x_ = -5e3:dx:5e3;
y_ = -5e3:dy:5e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

% Min = 0.0001;
% Max = 1;
% 
% R0 = Max .* (rand(size(x)) + 1j .* rand(size(x)));
% 
% % R0 = (rand(size(x))./sqrt(2) + 1j .* rand(size(x))./sqrt(2));
% 
% 
% figure; imagesc(x_, y_, (abs(R0))); colormap('jet')
% 
% shading flat; colorbar; axis equal tight;

%% GENERATE RANDOM FUNCTION 
SNR_db = 40;
SNR = 10^(SNR_db/10);

f = 0.045;

% rand('seed', 0);

nx1 = 1 + floor(0.9999.*f*length(x_));
nx2 = length(x_) - nx1 + 1; 

ny1 = 1 + floor(0.9999.*f*length(y_));
ny2 = length(y_) - ny1 + 1; 

F = (rand(length(x_), length(y_)) + 1j .* rand(length(x_), length(y_)) - 0.5);
F = fft2(F);
% figure; surface(x_.*1e-3, y_.*1e-3, db(abs(F).')./2); colorbar; shading interp;

F(nx1:nx2, :) = 0;
F(:, ny1:ny2) = 0;

F             = (ifft2(F));
N             = sqrt(max(abs(F)).^2./SNR) .* (rand(length(x_), length(y_)) + 1j .* rand(length(x_), length(y_)));
F             = (F + N)./max(max(F + N));

% pcolor(x_.*1e-3, y_.*1e-3, db(abs(F).')./2); colorbar; 

% pcolor(x_.*1e-3, y_.*1e-3, db(abs(N).')./2); colorbar; 
% shading interp; axis equal tight; 
% 
% xlim([-5 5]);
% ylim([-5 5]);
% 
% T = [-5:5:5];
% L = {};
% for m = 1:length(T)
%     if T(m) == 0
%         L{m} = '0';
%     else
%         L{m} = num2str(T(m), '%3.1f');
%     end
% end
% 
% set(gca, 'XTick', T, 'XTickLabel', L);
% set(gca, 'YTick', T, 'YTickLabel', L);
% 
% set(gca, 'FontSize', 18, 'LineWidth', 4);
% 
% title('Relectivity $\eta(x, y)$', 'Interpreter', 'latex', 'FontSize', 22);
% xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);

%% VELOCITY FIELDS


U = normrnd(10, 0.1, size(x));
% <<<<<<< Updated upstream
% V = normrnd(0, 0.1, size(x));
% =======
V = normrnd(10, 0.1, size(x));
% >>>>>>> Stashed changes
% Lx = x(end);
% Ly = y(end);
% 
% U = 2 + sin(2 * pi* x/Lx);
% V = 3 + sin(2 * pi* y/Ly);


PRT = 1;

t = eps:PRT:512*PRT;

Rt = zeros(length(t), length(y), length(x));
Rt(1, :, :) = F;

ETa_energy(1) = sum(sum(squeeze(abs(Rt(1, :, :)).^2) .* dx .* dy));
% 
% myMovie = zeros(1,length(t));

for k  = 1:length(t) - 1
    for i = 1:length(y)
        for l = 1:length(x)
           if ((i == 1) && (l == 1))
                Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l));
           elseif ((i == 1) && (l == length(x)))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l));
           elseif ((i == length(y)) && (l == 1))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l));
           elseif ((i == length(y)) && (l == length(x)))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l));
           
           
           elseif ((i == 1)) && (l ~= 1) && (l ~= length(x))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - V(i, l)./(2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l));
           elseif ((l == 1)) && ((i ~= 1) && (i ~= length(y)))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - V(i, l) ./( 2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l));
           
           elseif ((i == length(y))) && (l ~= 1) && (l ~= length(x))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - V(i, l)./(2 .* dy) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l));
           elseif ((l == length(x))) && (i ~= 1) && (i ~= length(y))
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - V(i, l)./( 2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l));
           else
               Rt(k + 1, i, l) = Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - V(i, l)./( 2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l));
           end
          
        end
    end
    
    ETa_energy(k+1) = sum(sum(squeeze(abs(Rt(k+1, :, :)).^2) .* dx .* dy));
%     surface(x_, y_, abs(squeeze(Rt(k+1, :, :)))); shading flat; colorbar; colormap('jet');
% %    colormap()
%     pause(0.1);
%     
end

close all;

Enable_movie = input('Do you want yo plot the movie? [0]: , ');
h = surface(x_, y_, abs(squeeze(Rt(1, :, :))));

if Enable_movie == 1
    for k = 2:length(t)
%         tic;
        
        h.CData = abs(squeeze(Rt(k, :, :)));
        
        shading flat; colormap('jet'); colorbar; % caxis([0.9 1]);
        pause(0.1);
        caption = sprintf('Frame #%d of %d, t = %.1f', k, length(t), t(k));
        title(caption, 'FontSize', 15); 
%         a = toc; 
%         disp(a);
        
    end
end

%movie(myMovie);


% figure; plot(t, ETa_energy);
