clear;
close all;

% Random Fields in 2D


dr = 227;
dph = pi/180;
r = 200:dr:15e3;
ph = eps:dph:2*pi;

[rho, phi] = meshgrid(r, ph); % Meshgrid of rho and phi forming a circle 

% GENERATE RANDOM FUNCTION 
SNR_db = 40;
SNR = 10^(SNR_db/10);

f = 0.02;

% rand('seed', 0);

nr1 = 1 + floor(0.9999.*f*length(r));
nr2 = length(r) - nr1 + 1; 

np1 = 1 + floor(0.9999.*f*length(ph));
np2 = length(ph) - np1 + 1; 



F = 10^4 .* (rand(length(r), length(ph)) + 1j .* rand(length(r), length(ph)) - 0.5);
figure; pcolor(rho .* cos(phi) .* 1e-3, rho .* sin(phi) .* 1e-3, db(abs(F).')./2); colorbar; shading interp;

F = fft2(F);
F(nr1:nr2, :) = 0;
F(:, np1:np2) = 0;

F             = (ifft2(F));
N             = sqrt(max(abs(F)).^2./SNR) .* (rand(length(r), length(ph)) + 1j .* rand(length(r), length(ph)));
F             = F + N;

figure; pcolor(rho .* cos(phi) .* 1e-3, rho .* sin(phi) .* 1e-3, db(abs(F).')./2); colorbar; shading interp;

% pcolor(x_.*1e-3, y_.*1e-3, db(abs(N).')./2); colorbar; 
shading interp; 

% axis equal tight; 

xlim([-15 15]);
ylim([-15 15]);

% T = [-15:5:15];
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

set(gca, 'FontSize', 18, 'LineWidth', 4);

title('Relectivity $\eta(x, y)$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);

%% 
%% R1 is the cloud reflectivity field
PRT = 100;

t = eps:PRT:8*PRT;

Rt = zeros(length(t), length(ph), length(r));
Rt(1, :, :) = F.';

U = normrnd(8, 1, size(rho));
V = normrnd(0, 0, size(rho));

Vr = U .* cos(phi) + V .* sin(phi);
Vaz = U .* sin(phi) - V .* cos(phi);

for k  = 1:length(t) - 1
    for i = 1:length(ph)
        for l = 1:length(r)
           if ((i == 1) && (l == 1))
                Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ (2 .* dr) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - Vaz(i, l) ./ r(l)./ (2 .* dph) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l));
           elseif ((i == 1) && (l == length(r)))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ (2 .* dr) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - Vaz(i, l) ./ r(l)./ (2 .* dph) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l));
           elseif ((i == length(ph)) && (l == 1))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ (2 .* dr) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - Vaz(i, l) ./ r(l)./ (2 .* dph) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l));
           elseif ((i == length(ph)) && (l == length(r)))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ (2 .* dr) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - Vaz(i, l) ./ r(l)./ (2 .* dph) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l));
           
           
           elseif ((i == 1)) && (l ~= 1) && (l ~= length(r))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ ( 2.* dr ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - Vaz(i, l) ./ r(l)./(2 .* dph) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l));
           elseif ((l == 1)) && ((i ~= 1) && (i ~= length(ph)))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ ( 2.* dr ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - Vaz(i, l) ./ r(l)./( 2 .* dph) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l));
           
           elseif ((i == length(ph))) && (l ~= 1) && (l ~= length(r))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ ( 2.* dr ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - Vaz(i, l) ./ r(l)./(2 .* dph) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l));
           elseif ((l == length(r))) && (i ~= 1) && (i ~= length(ph))
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ ( 2.* dr ) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - Vaz(i, l) ./ r(l)./( 2 .* dph) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l));
           else
               Rt(k + 1, i, l) = Rt(k, i, l) - Vr(i, l) ./ ( 2.* dr ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - Vaz(i, l) ./ r(l)./( 2 .* dph) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l));
           end
        end
    end
%     figure(k); surface(rho .* cos(phi), rho .* sin(phi), db(abs(squeeze(Rt(k+1, :, :))))/2); ...
%     shading flat; colormap('jet'); colorbar;
end




% figure;surface(rho .* cos(phi), rho .* sin(phi), db(abs(squeeze(Rt(3, :, :))))); ...
%     shading flat; colormap('jet'); colorbar;

makevideo(Rt, rho, phi, t)
