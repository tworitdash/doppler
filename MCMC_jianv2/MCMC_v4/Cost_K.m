function [E_samos] = Cost_K(k, U)

Udown = U(1:end-1, 1:k);
Uup = U(2:end, 1:k);



% Phi_hat = pinv(Udown) * Uup;
% 
% E = sum(sum(abs((Udown * Phi_hat - Uup))).^2);


Uupdown = [Uup Udown];

[Y G W] = svd(Uupdown);


end