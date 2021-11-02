function [a, b] = FourierCoeff(Vr, Nf, phi_axis)

    Nphi = length(phi_axis);
    a = zeros(1, Nf+1);
    b = zeros(1, Nf+1);
    for i = 1:Nf+1
        a(i) = 2/Nphi .* sum(Vr .* cos((i - 1) .* phi_axis));
        b(i) = 2/Nphi .* sum(Vr .* sin((i - 1) .* phi_axis));
    end

end