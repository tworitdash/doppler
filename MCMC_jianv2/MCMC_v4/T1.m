clear; 

% N = 10000; 
% Ns = 1000;
% 
% for i = 1:N
%     y(i, :) = rand(1, Ns);
% end
% 
% ysum = sum(y, 1);
% 
% histogram(ysum);


%% 

Mu = 0;
Sigma = 1;

X = linspace(-10, 10, 10000);

Y = 1/sqrt(2*pi*Sigma.^2) * exp(-(X - Mu).^2./(2*Sigma.^2));

Z = 1/(2) .* (erf((2*pi-X+Mu)./(sqrt(2))) - erf((-X+Mu)./(sqrt(2))));

figure; plot(X, Y); hold on; plot(X, Z);