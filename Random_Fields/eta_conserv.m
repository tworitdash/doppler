
function [J] = eta_conserv(eta, eta0, Vr, Vaz, rho, dr, dph, PRT)


       E = (eta - eta0)./PRT + Vr .* (eta - eta0)./dr + 1./rho .* Vaz .* (eta - eta0)./dph;
      
       J = sum(sum(abs(E).^2 .* rho .* dr .* dph));
end