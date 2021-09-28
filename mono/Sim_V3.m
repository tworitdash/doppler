%% Simulation of radar signals for weather application

K = 0.93; % Water droplet ratio of refractive index

dR = 2.27; 
R_max = 15e3;

r = eps:dR:R_max;

phi = linspace(eps, 360, 360) .* pi/180;

[R, Phi] = meshgrid(r, phi);

% figure; surface(r, phi, ones(size(R))); shading flat; colorbar;

X = normrnd(1, 0.5, size(R)) + normrnd(10, 1, size(R)); % + 19 + (20 -19) .* rand(size(R));

figure; surface(r, phi, X); shading flat; colorbar;

%% Random field generator


% corr.name = 'gauss';
% corr.c0 = 1;
% corr.c1 = 0; 
% corr.sigma = 2;
% mesh1 = R;
% 
% C = correlation_fun(corr,mesh1);


