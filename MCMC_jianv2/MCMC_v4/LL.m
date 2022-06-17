function LLout = LL(input_info, out)
%% Likelihood function of the radar data


z = out.Z_avail_vec;
t_avail = input_info.t_avail;

%% Geometrical information again evaluated 

% [geo] = Geo(input_info);
% 
xc0 = input_info.xc0;
yc0 = input_info.yc0;
rc0 = input_info.rc0;

% Nt = (input_info.N_gap + input_info.N_pulse) * input_info.N_rot; % number of time samples in HD signal
Navail = length(input_info.t_avail);

x(1) = sum( exp(1j .* (rc0) .* 4*pi/input_info.RADAR.lambda) );
   
if input_info.velocity.type == 1 % 1 for config with drops
%     D = linspace(input_info.velocity.D_min, input_info.velocity.D_max, input_info.NScatters);
% 
%     dD = D(2) - D(1);
%     ND = input_info.velocity.N0 * exp(-3.67 * D./input_info.velocity.D0);
%     D_int = ND.*D.^6.*dD;
%     Wt = 9.65 - 10.3 .* exp(-600 .* D .* 1e-3)
else
     u = normrnd(input_info.velocity.u.mu, input_info.velocity.u.sigma,[1 input_info.NScatters]);
     v = normrnd(input_info.velocity.v.mu, input_info.velocity.v.sigma, [1 input_info.NScatters]);
end

    xc(1, :) = xc0; 
    yc(1, :) = yc0;
    rc(1, :) = rc0;

for ti = 2:Navail
   xc(ti, :) = xc(ti - 1, :) + u .* (t_avail(ti) - t_avail(ti - 1));
   yc(ti, :) = yc(ti - 1, :) + v .* (t_avail(ti) - t_avail(ti - 1));
   rc(ti, :) = sqrt(xc(ti, :).^2 + yc(ti, :).^2);
    
   x(ti) = sum( exp(1j .* rc(ti, :) .* 4 .* pi ./ input_info.RADAR.lambda) );
end

xri = [real(x), imag(x)];
zri = [real(z), imag(z)];

% LLout = -Navail * log(2*pi) - Navail*2*log(out.sigma_n) -([zri - xri] * [zri - xri].')./(2 .* out.sigma_n^2);

LLout =  -([zri - xri] * [zri - xri].')./(2 .* out.sigma_n^2);




% end