function [Rt] = RFT(x_, y_, x, y, dx, dy, R1, U, V, PRT, Nt)

t = eps:PRT:Nt*PRT;

tau = t(end);
P = tau./(t + PRT);

Rt = zeros(length(t), length(y_), length(x_));
Rt(1, :, :) = R1;

ETa_energy(1) = sum(sum(squeeze(abs(Rt(1, :, :)).^2) .* dx .* dy));

for k  = 1:length(t) - 1
    for i = 1:length(y_)
        for l = 1:length(x_)
           if ((i == 1) && (l == 1))
                Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l)));
           elseif ((i == 1) && (l == length(x)))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l)));
           elseif ((i == length(y)) && (l == 1))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l)));
           elseif ((i == length(y)) && (l == length(x)))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ (2 .* dx) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - V(i, l)./ (2 .* dy) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l)));
           
           
           elseif ((i == 1)) && (l ~= 1) && (l ~= length(x))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - V(i, l)./(2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i, l)));
           elseif ((l == 1)) && ((i ~= 1) && (i ~= length(y)))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l))...
               - V(i, l) ./( 2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l)));
           
           elseif ((i == length(y))) && (l ~= 1) && (l ~= length(x))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - V(i, l)./(2 .* dy) .* PRT .* (Rt(k, i, l) - Rt(k, i - 1, l)));
           elseif ((l == length(x))) && (i ~= 1) && (i ~= length(y))
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l) - Rt(k, i, l - 1))...
               - V(i, l)./( 2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l)));
           else
               Rt(k + 1, i, l) = P(k) .* (Rt(k, i, l) - U(i, l) ./ ( 2.* dx ) .* PRT .* (Rt(k, i, l + 1) - Rt(k, i, l - 1))...
               - V(i, l)./( 2 .* dy) .* PRT .* (Rt(k, i + 1, l) - Rt(k, i - 1, l)));
           end
          
        end
    end
    
    ETa_energy(k+1) = sum(sum(squeeze(abs(Rt(k+1, :, :)).^2) .* dx .* dy));
    
%     figure; surface(x_, y_, db(abs(squeeze(Rt(k+1, :, :))))); ...
%     shading flat; colormap('jet'); colorbar;
end

end