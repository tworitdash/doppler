clear;
close all;

%% Random Fields in 2D


dx = 2270;
dy = 2270;
x_ = -15e3:dx:15e3;
y_ = -15e3:dy:15e3;

[x, y] = meshgrid(x_, y_); % Meshgrid of rho and phi forming a circle 

% Min = 0.0001;
Max = 1;

R0 = Max .* (rand(size(x)) + 1j .* rand(size(x)));

% R0 = (rand(size(x))./sqrt(2) + 1j .* rand(size(x))./sqrt(2));


figure; imagesc(x_, y_, (abs(R0))); shading flat;colorbar; axis equal tight;

lambda = 1e3;

k0 = 2*pi/lambda;

Kx_ = linspace(-k0, k0, 100);
dKx = Kx_(2) - Kx_(1);

Ky_ = linspace(-k0, k0, 100);
dKy = Ky_(2) - Ky_(1);

[Kx, Ky] = meshgrid(Kx_, Ky_);


for i = 1:length(Kx_)
    for l = 1:length(Ky_)

        Spec(l, i) = sum(sum(R0 .* exp(1j .* (Kx_(i) .* x + Ky_(l) .* y)) .* dKx .* dKy));
    
    end
end


figure; surface(Kx, Ky, abs(Spec)); shading flat; colormap('jet'); colorbar;


Spec_modified = Spec .* log .* exp(1j .* 2 * pi * randn(size(Kx)));

% figure; surface(K .* cos(Th), K .* sin(Th), abs(Spec_modified)); shading flat; colormap('jet'); colorbar;

for i = 1:length(x_)
    for l = 1:length(y_)
        R1(l, i) = sum(sum(Spec_modified .* exp(1j .* (Kx .* x_(i) + Ky .* y_(l))) .* dx .* dy));
    end
end


figure; surface(x, y, db(abs(R1))/2); shading flat; colormap('jet'); colorbar;


for i = 1:length(Kx_)
    for l = 1:length(Ky_)

        Spec_R1(l, i) = sum(sum(R1 .* exp(1j .* (Kx_(i) .* x + Ky_(l) .* y)) .* dKx .* dKy));
    
    end
end
figure; surface(Kx, Ky, abs(Spec_R1)); shading flat; colormap('jet'); colorbar;




%% R1 is the cloud reflectivity field
PRT = 10;

t = eps:PRT:100*PRT;

Rt = zeros(length(t), length(phi), length(r));
Rt(1, :, :) = R1;

U = 8;
V = 8;

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
           if db(abs(Rt(k+1, i, l))) > 40
               Rt(k + 1, i, l) = 10^(40/20);
           end
        end
    end
%     figure(k); surface(rho .* cos(phi), rho .* sin(phi), db(abs(squeeze(Rt(k+1, :, :))))); ...
%     shading flat; colormap('jet'); colorbar;
end




% figure;surface(rho .* cos(phi), rho .* sin(phi), db(abs(squeeze(Rt(3, :, :))))); ...
%     shading flat; colormap('jet'); colorbar;

makevideo(Rt, rho, phi, t)