clear;
close all;

%% Gaussian Random Fields in 2D
dr = 227;
dph = 1*pi/180;
r = 100:dr:15e3;
ph = eps:dph:2*pi;

[rho, phi] = meshgrid(r, ph); % Meshgrid of rho and phi forming a circle 


R0 = rand(size(rho)) + 1j .* rand(size(rho)) - 0.5;

% figure; surface(rho .* cos(phi), rho .* sin(phi), abs(R0)); shading flat; colormap('jet'); colorbar;

lambda = 2e3;

k0 = 2*pi/lambda;

K_rho = linspace(eps, k0, 100);
dK = K_rho(2) - K_rho(1);

Theta = linspace(eps, 2*pi, 360);
dTheta = Theta(2) - Theta(1);

[K, Th] = meshgrid(K_rho, Theta);


for i = 1:length(K_rho)
    for l = 1:length(Theta)

        Spec(l, i) = sum(sum(R0 .* exp(1j .* K_rho(i) .* rho .* cos(phi - Theta(l))) .* rho .* dr * dph));
    
    end
end

figure; surface(K .* cos(Th), K .* sin(Th), abs(Spec)); shading flat; colormap('jet'); colorbar;


% Spec_modified = Spec .* exp(K.^2 .* cos(Th) +K .* sin(Th)+3) .* ((cos(Th)).^2) .* exp(1j .* 2 * pi * randn(size(K)));

f = 0.03;

nr1 = 1 + floor(0.9999.*f*length(r));
nr2 = length(r) - nr1 + 1; 

np1 = 1 + floor(0.9999.*f*length(ph));
np2 = length(ph) - np1 + 1; 


Spec_modified = Spec;
Spec_modified(np1:np2, :) = 0;
Spec_modified(:, nr1:nr2) = 0;
 
% figure; surface(K .* cos(Th), K .* sin(Th), abs(Spec_modified)); shading flat; colormap('jet'); colorbar;

for i = 1:length(r)
    for l = 1:length(ph)
        R1(l, i) = sum(sum(Spec_modified .* exp(-1j .* K .* r(i) .* cos(ph(l) - Th)) .* K .* dK * dTheta));
    end
end


figure; surface(rho .* cos(phi), rho .* sin(phi), db(abs(R1))); shading flat; colormap('jet'); colorbar;


for i = 1:length(K_rho)
    for l = 1:length(Theta)

        Spec_R1(l, i) = sum(sum(R1 .* exp(1j .* K_rho(i) .* rho .* cos(phi - Theta(l))) .* rho .* dr * dph));
    
    end
end
% figure; surface(K .* cos(Th), K .* sin(Th), abs(Spec_R1)); shading flat; colormap('jet'); colorbar;




%% R1 is the cloud reflectivity field
PRT = 10;

t = eps:PRT:8*PRT;

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