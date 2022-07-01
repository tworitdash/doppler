clear; 
close all;
Ns = 100000;
lambda = 3e-2;
dT = 1e-3;

sig = 1 * 4 * pi/lambda * dT;
mu = 7 * 4 * pi/lambda * dT;


u = mod(normrnd(mu, sig, [1 Ns]), 2 * pi);
k = 10;
figure; histogram(u*lambda/(4 * pi * dT));


% v = normrnd(cos(k * mu), k .* sin(k * mu).*sig, [1 Ns]);


figure; histogram((cos(mod(k * u, 2 * pi)))); 

%hold on; histogram(v);

y = linspace(-1, 1, 10000);

Fy = 1./(sqrt(1 - y.^2)) .* 1./sqrt(2 * pi * sig.^2)  .* exp(-(acos(y) - mu).^2./(2 * sig.^2));

figure; plot(y, Fy);


%% 

mu = pi/4;
sigma = pi/10;

x = normrnd(mu, sigma, [1 1000000]);

n = normpdf(x, mu, sigma);

Exp = sum(x .* n)/sum(n);

figure; plot(n);

y = cos(x);

Expy = sum(y .* n)./sum(n);
Expyan = exp(-(sigma)^2/2).* cos(mu);

% Expyanm = sigma .* exp(-1j .* mu - sigma^2/4) ./ (2 * sqrt(2) .* sigma) .* ( -erfz(mu/sigma - 1j .* sigma/2 - pi/sigma) ...
%     + erfz(mu/sigma - 1j .* sigma/2 + pi/sigma) + exp(1j .* 2 .* mu) ...
%     .* ( -erfz(mu/sigma + 1j .* sigma/2 - pi/sigma) + erfz(mu/sigma + 1j .* sigma/2 + pi/sigma) ) );



Expyam = 1./(2 .* sigma .* sqrt(sigma^2)) .* exp(-1j .* mu - 1/(2 .* sigma^2)) .* ...
    ( -erfz(mu/sigma - (1j .*  sigma)/2 - pi/sigma) + ... 
   erfz(mu/sigma - (1j .* sigma)/2 + pi/sigma) + ...
   exp(2 .* 1j .* mu) .* (-erfz(mu/sigma + 1j .* sigma/2 - pi/sigma) + ...
      erfz(mu/sigma + (1j .* sigma)/2 + pi/sigma))) .* ...
      ...
  (erfz((1j + sigma^2 .* (-mu + pi))/(sqrt(2).* sigma)) - erfz((1j - sigma^2 .* (mu + pi))/(sqrt(2) .* sigma)) - ... 
   1j .* exp(2 .* 1j .* mu) .* (erfi((1 + 1j .* sigma^2 .* (-mu + pi))/(sqrt(2) .* sigma)) - ...
      erfi((1 - 1j .* sigma^2 .* (mu + pi))./(sqrt(2) .* sigma)))); 