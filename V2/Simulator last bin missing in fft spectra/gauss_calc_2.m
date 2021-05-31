function [M0t,M1t,M2t] = gauss_calc_2(x,v,dv)
% dv = 2*vmax/(Ndoppler);
% v = -vmax:dv:vmax-dv;

M0t = sum(x*dv);

if M0t == 0
    M1t = NaN;
    M2t = NaN;
else
M1t = (1/M0t)*sum(v.*abs(x).*dv);
M2t = sqrt((1/(M0t)*sum((v-M1t).^2.*abs(x).*dv)));
end