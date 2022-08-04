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


Nt = 128;
K = 0:1:Nt-1;
% K = Nt

x = linspace(-2*pi*Nt, 2*pi*Nt, 10000);
dx = x(2) - x(1);

for k = 1:length(K)
    
%      if mod(K(k), 10) == 0
%       figure(1000); hold on; histogram(K(k) .* U + beta);
%      end
     
     sig_(k, :) = ( exp(1j .* K(k) .* U + 1j .* beta) );
     sig_(k, :) = sig_(k, :);
     sig(k) = sum(sig_(k, :));
%      sig(k) = sum( exp(1j .* K(k) .* U ) );
%     
    
end

 % hold on; histogram(imag(sig_(1, :)));

%% CDF Technique
% s = rand(1, Nu);
% 
% for  k = 1:length(K)
%     F(k, :) = 1/(4 * pi) .* ( (K(k).*Mu - acos(s)) .* erf( (K(k) .* Mu - acos(s))./(sqrt(2) .* K(k) .* Sigma) ) + ...
%         (-K(k).*Mu + acos(s) - 2*pi) .* erf( (K(k) .* Mu - acos(s) + 2*pi)./(sqrt(2) .* K(k) .* Sigma) )  + ...
%         sqrt(2./pi) .* K(k) .* Sigma .* ( exp(-(acos(s) - K(k) .* Mu).^2./(2.*K(k).^2.*Sigma.^2))...
%         - exp(-(acos(s) - K(k).*Mu + 2*pi).^2./(2.*K(k).^2.*Sigma.^2)) ) );
% end
% 
% 
% 
% figure; histogram((F(1, :))); %hold on; histogram(imag(F(10, :)));

s = linspace(-1+1e-5, 1-1e-5, Nu);% s = rand(1, Nu);

for k = 1:length(K)

    f(k, :) = -1./(4 * pi * sqrt(1 - s.^2)) .* ( erf( (K(k).*Mu - acos(s))./(sqrt(2).*K(k).*Sigma) ) - ...
        erf( (K(k).*Mu - acos(s) + 2 * pi)./(sqrt(2).*K(k).*Sigma) ) );% .* sqrt(Nu); 
    
end

figure; histogram(real(sig_(12, :))); hold on; plot(s, real(f(12, :))); % hold on; histogram(imag(sig_(1, :)));