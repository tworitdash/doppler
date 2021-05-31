function f = gauss_gen_2(v,M0,M1,M2)
% dv = 2*vmax/(Ndoppler);% velocity resolution
% x = -vmax:dv:vmax-dv;% velocity bins

f = M0/(sqrt(2*pi)*M2)*exp(-(v-M1).^2/(2*M2^2));