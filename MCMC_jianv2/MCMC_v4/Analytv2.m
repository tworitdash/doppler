
clear;
% close all;

lambda = 0.03;
dt = 1e-3;

Mu = 7.5 .* 0.4189;
Sigma = 2 .*  0.4189;

% rng(1, 'twister');
Nu = 100000;

beta = 2 * pi .* rand([1 Nu]); 

    
U = normrnd(Mu, Sigma, [1 Nu]);


Nt = 1024;
K = 0:1:Nt-1;
% K = Nt

x = linspace(-2*pi*Nt, 2*pi*Nt, 10000);
dx = x(2) - x(1);

for k = 1:length(K)
    
     if mod(K(k), 10) == 0
      figure(1000); hold on; histogram(K(k) .* U + beta);
     end
     
     sig(k) = sum( exp(1j .* K(k) .* U + 1j .* beta) );
%      sig(k) = sum( exp(1j .* K(k) .* U ) );
%     
    
end

figure; plot(real(sig)); hold on; plot(imag(sig));

figure(100); hold on; histogram(real(sig)); figure(200); hold on; histogram(imag(sig));

Mur = mean(real(sig)); Mui = mean(imag(sig));
Sigr = std(real(sig)); Sigi = std(imag(sig));

% 


% Mucos = 1./(8 .* pi) .* exp(-(1/2 .* K .* (2 .*1j .* Mu + K .* Sigma^2))) .* (-1j .* exp(2 .* 1j .* k .* Mu) .* ...
%      erfz((k .* Mu + pi + 1j .* K.^2 .* Sigma^2)./(sqrt(2) .* K .* Sigma)) - ...
%    1j .* erfz((pi + K .* (Mu - 1j .* K .* Sigma^2))./(sqrt(2) .* K .* Sigma)) + ...
%    1j .* exp(2 .* 1j .* K .* Mu) .* erfz((pi + K .* (Mu + 1j .* K .* Sigma.^2))./(sqrt(2) .* K .* Sigma)) + ...
%    exp(2 .* 1j .* K .* Mu) .* erfi((-1j .* Mu + K .* Sigma^2)./(sqrt(2) .* Sigma)) + ...
%    erfi( (1j .* Mu + K .* Sigma^2)/(sqrt(2) .* Sigma) ) + ...
%    erfi( (1j .* K .* Mu + 1j .* pi + K.^2 .* Sigma^2)./(sqrt(2) .* K .* Sigma) ) - ... 
%    exp(2 .* 1j .* K .* Mu) .* erfi((-1j .* K .* Mu - 2 .* 1j .* pi + K.^2 .*  Sigma^2)./(sqrt(2) .* K .* Sigma)) - ...
%    erfi((1j .* K .* Mu + 2 .* 1j .* pi + K.^2 .* Sigma^2)./(sqrt(2) .* K .* Sigma)));

