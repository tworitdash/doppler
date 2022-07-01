%%
clear;

lambda = 3e-2;
dt = 1e-3;

mu = 5*4*pi/lambda*dt;
sigma = 0.5 * 4*pi/lambda*dt;
Nu = 100000;

u = linspace(0, 15, Nu)*4*pi/lambda*dt;

K = [0:1:4 101:1:105 201:1:205];

Nt = length(K);


for i = 1:Nu

    F(i) = -sum(1 + exp(-K.^2*sigma^2) - 2 .* exp(-K.^2.*sigma^2/2) .* cos(K .* (mu - u(i))));
    
end

figure; plot(u/(4*pi/lambda*dt), F); grid on;