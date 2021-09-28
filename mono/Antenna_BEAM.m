function [wtx] = Antenna_BEAM(D, lambda, ph, ph_0)
    
    wtx_num = 8 .* (besselj(2, pi * D .* sin(ph - ph_0)/lambda)).^2;
    wtx_denum = (pi * D .* sin(ph - ph_0)/lambda).^4;
    
    wtx = wtx_num./wtx_denum;
end