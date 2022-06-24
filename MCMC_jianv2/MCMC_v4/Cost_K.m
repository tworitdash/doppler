function [E_samos] = Cost_K(k, U)

Udown = U(1:end-1, 1:k);
Uup = U(2:end, 1:k);



% Phi_hat = pinv(Udown) * Uup;
% 
% E = sum(sum(abs((Udown * Phi_hat - Uup))).^2);


Uupdown = [Uup Udown];

[~, G, ~] = svd(Uupdown);

gamma = diag(G);

E_samos = 1/k * sum(gamma(k+1:end));


end