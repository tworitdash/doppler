function dJmdZ = besselj_der(m, z)

    %dz = 1e-9;
    %dJmdZ = (besselj(m, z + dz) - besselj(m, z - dz))./(2 * dz);
    
    
    dJmdZ= (besselj(m - 1, z) - besselj(m + 1, z))./2;
end