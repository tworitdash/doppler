clear;
lambda = 3e-2;
dT = 1e-3;
Nt = 5;

t = eps:dT:Nt*dT;

Nu = 100000;

u_axis = linspace(0, 15, Nu) * (4 * pi)/lambda * dT;

Ns = 10000;

Mu_u = 5 * (4 * pi)/lambda * dT;
Sigma_u = 1 * 4 * pi / lambda * dT;

U = normrnd(Mu_u, Sigma_u, [1 Ns]);


% 
% for m = 1:Nu
% 
% for k = 1:Nt
% 
% for i = 1:Ns
%     
%     ui = U(i);
% 
%     y = @(um) cos((ui - um) .* k);
% 
%     Muy = y(Mu_u);
%     Sigyy = sqrt(k.^2 .* Sigma_u^2 .* sin((ui - Muy).*k)^2); 
% 
%     Int1 = @(um) y(um) .* exp(-(y(um) - Muy).^2./(2 .* Sigyy.^2)) ./sqrt(2 .* pi .* Sigyy.^2);
% 
%     Int1_val = Int1(U);
% 
%     Int1_sum(i) = sum(Int1_val);
% 
% end
% 
% 
% 
% First_Term(k) = sum(Int1_sum);
% 
% 
% 
% Tau = @(um) cos(um - (u_axis(m)) .* k);
% 
% 
% %% 
% 
% 
% 
% Mtau = Tau(Mu_u);
% Sigtau = sqrt(k.^2 .* Sigma_u^2 .* sin((Muy - u_axis(m)).*k)^2); 
% 
% Int2 = @(um) Tau(um) .* exp(-(Tau(um) - Mtau).^2./(2 .* Sigtau.^2)) ./sqrt(2 .* pi .* Sigtau.^2);
% Int2_val = Int2(U);
% 
% Second_Term(k) = sum(Int2_val);
% 
% end
% 
% LL(m) = Nt .* (Ns + 1) + 2 .* sum(First_Term) - 2 .* sum(Second_Term);
% 
% end
% 
% 
% figure; plot(u_axis*(lambda)/(4 * pi * dT), LL);




%% 


for m = 1:Nu
    
    for k = 1:Nt
        LLk(k) = -( (sum(cos(U .* k)) - cos(u_axis(m).*k)).^2 + ...
            (sum(sin(U .* k)) - sin(u_axis(m).*k)).^2 );
    end
    
    LL(m) = sum(LLk);
end

figure; plot(u_axis*(lambda)/(4 * pi * dT), LL);