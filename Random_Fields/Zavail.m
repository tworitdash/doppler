function [zi, z] = Zavail(z_model, n_rot, N_Sweep, N_sec)

z = reshape(z_model, [N_sec * N_Sweep, n_rot]);

for i = 1:n_rot
    zi(i, :) = z(1:N_Sweep, i);
end

% figure; plot(real(z_model)); hold on; plot(real(zi(1, :))); 
% figure; plot(imag(z_model)); hold on; plot(imag(zi(1, :)));

end