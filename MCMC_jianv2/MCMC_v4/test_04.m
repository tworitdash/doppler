

Np = 5;
M = 2;
dT = 1e-3;
lambda = 0.03;

Tscan = 100 * dT;

ut = 5 * 4 * pi / lambda;

u = linspace(0, 15, 10000) * 4 * pi / lambda + 1e-7;

Y1 = sin(Np * dT * (ut - u)/2) ./ sin((ut - u)/2 * dT);
Y2 = sin(M .* Tscan .* (ut - u)/2) ./ sin(Tscan .* (ut - u)/2);
Y3 = cos((ut - u)/2 .* ((Np - 1) * dT + Tscan * (M - 1)));

% figure; plot(u, Y1./max(Y1)); hold on; plot(u, Y2./max(Y2)); 
% 
% figure; hold on; plot(u, Y1.*Y2./max(Y1.*Y2)); hold on; plot(u, Y1.* Y2 .* Y3 ./ max(Y1 .* Y2 .* Y3));

% figure; plot(u, Y1./max(Y1)); 



Y = 2 .* Y1 .* Y2 .* Y3 - 2 * M * Np;
hold on; plot(u*lambda/(4*pi), Y); grid on;
 
% X = sin(Np * (ut - u)/2) ./ sin((ut - u)/2) .* sin(Nscan .* Tscan .* (ut - u)/2) ./ sin(Tscan .* (ut - u)/2)...
%     .* cos((ut - u)/2 .* (Np + Tscan * (M - 1)));
% 
% figure; plot(u, Y./max(Y)); hold on; plot(u, (X./max(X))); grid on;


%% V2


%%
% dT = 1e-3;
% Np = 25*dT;
% M = 3;
% 
% Tscan = 100*dT;
% lambda = 0.03;
% 
% ut = 3 * 4 * pi/lambda;
% 
% u = ut + linspace(-ut, 15*4*pi/lambda-ut, 1000) +  eps;
% % u(u==ut) = ut+eps;
% 
% Y1 = sin(Np * (ut - u)/2 + eps) ./ sin((ut - u)/2 + eps);
% 
% Y2 = sin(M .* Tscan .* (ut - u)/2) ./ sin(Tscan .* (ut - u)/2);
% Y3 = cos((ut - u)/2 .* (Np + 1 + Tscan * (M - 1)));
% 
% % figure; plot(u, Y1./max(Y1)); hold on; plot(u, Y2./max(Y2)); 
% % 
% % figure; hold on; plot(u, Y1.*Y2./max(Y1.*Y2)); hold on; plot(u, Y1.* Y2 .* Y3 ./ max(Y1 .* Y2 .* Y3));
% 
% % figure; plot(u, Y1./max(Y1)); 
% hold on; plot(u*lambda/(4*pi), Y1.* Y2 .* Y3./max(Y1 .* Y2 .* Y3));
% 
% 
