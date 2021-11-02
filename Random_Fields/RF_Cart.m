clear;
close all;

%% Random Fields in 2D


dx = 100;
dy = 100;
x_ = eps:dx:15e3;
y_ = eps:dy:15e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

% Min = 0.0001;
Max = 1;

R0 = Max .* (rand(size(x)) + 1j .* rand(size(x)) - 0.5);

% R0 = (rand(size(x))./sqrt(2) + 1j .* rand(size(x))./sqrt(2));


% figure; imagesc(x_, y_, (abs(R0))); shading flat;colorbar; axis equal tight;

lambda = 1e3;

k0 = 2*pi/lambda;

Kx_ = linspace(-k0, k0, 100);
dKx = Kx_(2) - Kx_(1);

Ky_ = linspace(-k0, k0, 100);
dKy = Ky_(2) - Ky_(1);

[Kx, Ky] = meshgrid(Kx_, Ky_);


for i = 1:length(Kx_)
    for l = 1:length(Ky_)

        Spec(l, i) = sum(sum(R0 .* exp(-1j .* (Kx_(i) .* x + Ky_(l) .* y)) .* dKx .* dKy));
    
    end
end


% figure; surface(Kx, Ky, abs(Spec)); shading flat; colormap('jet'); colorbar;

f = 0.2;
Spec_modified = Spec;

nx1 = 1 + floor(0.9999.*f*length(Kx_));
nx2 = length(Kx_) - nx1 + 1; 

ny1 = 1 + floor(0.9999.*f*length(Ky_));
ny2 = length(Ky_) - ny1 + 1;


Spec_modified(nx1:end, :) = 0;
Spec_modified(:, ny1:end) = 0;



% figure; surface(Kx,Ky, abs(Spec_modified)); shading flat; colormap('jet'); colorbar;

for i = 1:length(x_)
    for l = 1:length(y_)
        R1(l, i) = sum(sum(Spec_modified .* exp(1j .* (Kx .* x_(i) + Ky .* y_(l))) .* dx .* dy));
    end
end


figure; surface(x, y, (abs(R1))/2); shading flat; colormap('jet'); colorbar;


for i = 1:length(Kx_)
    for l = 1:length(Ky_)

        Spec_R1(l, i) = sum(sum(R1 .* exp(-1j .* (Kx_(i) .* x + Ky_(l) .* y)) .* dKx .* dKy));
    
    end
end
% figure; surface(Kx, Ky, abs(Spec_R1)); shading flat; colormap('jet'); colorbar;




%% R1 is the cloud reflectivity field
U = normrnd(3, 0.1, size(x));
V = normrnd(1, 0.1, size(x));


PRT = 1e-3;

t = eps:PRT:100*PRT;

Rt = zeros(length(t), length(y_), length(x_));
Rt(1, :, :) = R1;

ETa_energy(1) = sum(sum(squeeze(abs(Rt(1, :, :)).^2) .* dx .* dy));

for k  = 1:length(t) - 1
    for i = 1:length(y_)
        for l = 1:length(x_)
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
    
%     figure; surface(x_, y_, db(abs(squeeze(Rt(k+1, :, :))))); ...
%     shading flat; colormap('jet'); colorbar;
end




% figure;surface(rho .* cos(phi), rho .* sin(phi), db(abs(squeeze(Rt(3, :, :))))); ...
%     shading flat; colormap('jet'); colorbar;

makevideo_XY(Rt, x_, y_, t)


%% RADAR Filter

lambda = 0.03;
R = 5e3;
dr = 227; dph = pi/180; 

phi = eps:dph:pi/2;
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
        F = squeeze(Rt(l, :, :));
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
            
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;
            
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);
            
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
            
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;
            
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);
            
            
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
            
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;
            
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);
            
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
            
            Rr(i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;
        
            A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));
            
            D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);
        
        end
    end
end

figure; surface(xr, yr, db(abs(Rr).')./2); colorbar; colormap('jet'); shading flat;

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
% figure; surface(xr, yr, Vr.'); colorbar; colormap('jet'); shading flat;

