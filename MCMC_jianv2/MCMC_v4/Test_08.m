clear;

lambda = 3e-2;
dt = 1e-3;

mu = 3*4*pi/lambda*dt;
sig = 1*4*pi/lambda*dt;
Nu = 100;
Nuobs = 10000;

u_ = linspace(0, 15, Nuobs)*4*pi/lambda*dt+1e-7;


U_ = normrnd(mu, sig, [1 Nu]); %*4*pi/lambda*dt+1e-7;

[u, U] = meshgrid(u_, U_);

[U1, U2] = meshgrid(U_, U_);

% [umn1, umn2] = meshgrid(U, U);

% K = [0:1:4 101:1:105 201:1:205 301:1:305];
K = [0:1:4];

% [u, U, K] = ndgrid(u_, U_, K_);


Nt = length(K);

Ns = 10000;

Np = 5;
M = 1;
Nscan = 100;

X = @(x, y) sin(Np .* (x - y)/2)./sin((x - y)./2) .* sin(Nscan .* M .* (x- y)/2)./sin(Nscan .* (x - y)./2) .* cos((x - y)./2.*(Np-1 + Nscan.*(M - 1)));


% F = 2 .* sin(Np .* (U - u)/2)./sin((U - u)./2) .* sin(Nscan .* M .* (U - u)/2)./sin(Nscan .* (U - u)./2) .* cos((U - u)./2.*(Np-1)) - 2 * Np .* M;

F = sum(sum(X(U1, U2))) + sum(X(U, u), 1);

figure(102); hold on; plot(u_./(4*pi/lambda*dt), F-max(F));


