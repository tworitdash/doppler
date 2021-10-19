clear;
close all;

%% Random Fields in 2D


dx = 10;
dy = 10;
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
SNR_db = 20;
SNR = 10^(SNR_db/10);


F = 10^3 .* (rand(length(x_), length(y_)) + 1j .* rand(length(x_), length(y_)) - 0.5 - 1j .* 0.5);
db(abs(min(min(F))))


figure; pcolor(x_.*1e-3, y_.*1e-3, (abs(F).')); colorbar; shading interp; 
shading interp; axis equal tight; 

xlim([-5 5]);
ylim([-5 5]);

set(gca, 'FontSize', 18, 'LineWidth', 4);

title('Relectivity $\eta(x, y)$ [Linear Scale]', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);

F = fft2(F);
% figure; surface(x_.*1e-3, y_.*1e-3, (abs(F).')./2); colorbar; shading interp;

f = 0.004;

nx1 = 1 + floor(0.9999.*f*length(x_));
nx2 = length(x_) - nx1 + 1; 

ny1 = 1 + floor(0.9999.*f*length(y_));
ny2 = length(y_) - ny1 + 1;


F(nx1:nx2, :) = 0;
F(:, ny1:ny2) = 0;
% 
F             = (ifft2(F));
figure; pcolor(x_.*1e-3, y_.*1e-3, (abs(F).')); colorbar; 
shading interp; axis equal tight; 

xlim([-5 5]);
ylim([-5 5]);


set(gca, 'FontSize', 18, 'LineWidth', 4);

title('Relectivity $\eta(x, y)$ [Linear Scale]', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);


% 
F             = imgaussfilt(real(F), 10);
N             = sqrt(max(max(abs(F))).^2./SNR) .* (rand(length(x), length(y)) + 1j .* rand(length(x), length(y)));
F             = F + N;

db(abs(min(min(F))))
figure; pcolor(x_.*1e-3, y_.*1e-3, (abs(F).')); colorbar; 

% pcolor(x_.*1e-3, y_.*1e-3, db(abs(N).')./2); colorbar; 
shading interp; axis equal tight; 

xlim([-5 5]);
ylim([-5 5]);


set(gca, 'FontSize', 18, 'LineWidth', 4);

title('Relectivity $\eta(x, y)$ [Linear Scale]', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);

%% RADAR Filter
lambda = 0.03;
R = 5e3;
dr = 227; dph = pi/180; 

phi = eps:dph:2 * pi;
r = eps:dr:R;


[rr, phir] = meshgrid(r, phi);

xr = rr .* cos(phir);
yr = rr .* sin(phir);

PRT = 1e-3; 
v_amb = lambda./(4 .* PRT);

u = 3; 
v = 1;

for i = 1:length(r)
    for l = 1:length(phi)
        
        if (i == length(r)) && (l == length(phi))
            xil1(i, l) = r(i) .* cos(phi(l));
            yil1(i, l) = r(i) .* sin(phi(l));
            xil2(i, l) = r(i) .* cos(phi(l));
            yil2(i, l) = r(i) .* sin(phi(l));

            [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
            [~, idy1(i, l)] = min(abs(x_ - yil1(i, l)));
            [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
            [~, idy2(i, l)] = min(abs(x_ - yil2(i, l)));
        
            idxmin(i, l) = min(idx1(i, l), idx2(i, l));
            idxmax(i, l) = max(idx1(i, l), idx2(i, l));
            
            idymin(i, l) = min(idy1(i, l), idy2(i, l));
            idymax(i, l) = max(idy1(i, l), idy2(i, l));
            Nxy(i, l) = (idxmax(i, l)-idxmin(i, l)+1) * (idymax(i, l)-idymin(i, l)+1);
            Axy(i, l) = Nxy(i, l) * dx .* dy;
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy./Axy(i, l);
            
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            delX = u * l * PRT;
            delY = v * l * PRT;
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delX).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delY).^2);
            
        elseif (i == length(r)) && (l ~= length(phi))
            xil1(i, l) = r(i) .* cos(phi(l));
            yil1(i, l) = r(i) .* sin(phi(l));
            xil2(i, l) = r(i) .* cos(phi(l + 1));
            yil2(i, l) = r(i) .* sin(phi(l + 1));

            [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
            [~, idy1(i, l)] = min(abs(x_ - yil1(i, l)));
            [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
            [~, idy2(i, l)] = min(abs(x_ - yil2(i, l)));

            idxmin(i, l) = min(idx1(i, l), idx2(i, l));
            idxmax(i, l) = max(idx1(i, l), idx2(i, l));
            
            idymin(i, l) = min(idy1(i, l), idy2(i, l));
            idymax(i, l) = max(idy1(i, l), idy2(i, l));
            Nxy(i, l) = (idxmax(i, l)-idxmin(i, l)+1) * (idymax(i, l)-idymin(i, l)+1);
            Axy(i, l) = Nxy(i, l) * dx .* dy;
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy./Axy(i, l);
            
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            delX = u * l * PRT;
            delY = v * l * PRT;
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delX).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delY).^2);
            
            
        elseif (i ~= length(r)) && (l == length(phi))
            xil1(i, l) = r(i) .* cos(phi(l));
            yil1(i, l) = r(i) .* sin(phi(l));
            xil2(i, l) = r(i + 1) .* cos(phi(l));
            yil2(i, l) = r(i + 1) .* sin(phi(l));

            [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
            [~, idy1(i, l)] = min(abs(x_ - yil1(i, l)));
            [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
            [~, idy2(i, l)] = min(abs(x_ - yil2(i, l)));

            idxmin(i, l) = min(idx1(i, l), idx2(i, l));
            idxmax(i, l) = max(idx1(i, l), idx2(i, l));
            
            idymin(i, l) = min(idy1(i, l), idy2(i, l));
            idymax(i, l) = max(idy1(i, l), idy2(i, l));
            
            Nxy(i, l) = (idxmax(i, l)-idxmin(i, l)+1) * (idymax(i, l)-idymin(i, l)+1);
            Axy(i, l) = Nxy(i, l) * dx .* dy;
            
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy ./ Axy(i, l);
            
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            delX = u * l * PRT;
            delY = v * l * PRT;
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delX).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delY).^2);
            
        else
            
            xil1(i, l) = r(i) .* cos(phi(l));
            yil1(i, l) = r(i) .* sin(phi(l));
            xil2(i, l) = r(i + 1) .* cos(phi(l + 1));
            yil2(i, l) = r(i + 1) .* sin(phi(l + 1));

            [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
            [~, idy1(i, l)] = min(abs(x_ - yil1(i, l)));
            [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
            [~, idy2(i, l)] = min(abs(x_ - yil2(i, l)));
            
            idxmin(i, l) = min(idx1(i, l), idx2(i, l));
            idxmax(i, l) = max(idx1(i, l), idx2(i, l));
            
            idymin(i, l) = min(idy1(i, l), idy2(i, l));
            idymax(i, l) = max(idy1(i, l), idy2(i, l));
            Nxy(i, l) = (idxmax(i, l)-idxmin(i, l)+1) * (idymax(i, l)-idymin(i, l)+1);
            Axy(i, l) = Nxy(i, l) * dx .* dy;
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy ./ Axy(i, l);
        
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            delX = u * l * PRT;
            delY = v * l * PRT;
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delX).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l)) + delY).^2);
        
        end
    end
end

figure; pcolor(xr.*1e-3, yr.*1e-3, db(abs(Rr).')./2); colorbar; shading flat; grid on;
% 
xlim([-5 5]);
ylim([-5 5]);


set(gca, 'FontSize', 18, 'LineWidth', 4);

title('Relectivity $Z(\rho, \phi) [dB]$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$X$ [km]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Y$ [km]', 'Interpreter', 'latex', 'FontSize', 22);

% for i = 1:length(r)
%     for l = 1:length(phi)
%         X = (A(i, l).a(:, :)).' .* exp(-1j .* 4 * pi ./ lambda .* D(i, l).d(:, :));
%         X_re = reshape(X, 1, [size(X, 1) .* size(X, 2)]);
%         XF = 1./sqrt(length(X_re)) .* fftshift(fft(X_re));
%         dv = lambda./(2 .* length(XF) .* PRT);
%         vel_axis = linspace(-v_amb, v_amb-dv, length(X_re));
%         PT(i, l) = sum(abs(XF).^2 .* dv);
%         Vr(i, l) = sum(vel_axis .* abs(XF).^2 .* dv)./PT(i, l);
%     end
% end
% 
% figure; surface(xr, yr, Vr.'); colorbar;  shading interp;
