function [PT, Mu, Sigma] = Ret(Rr, r, phi, N_sweep, PRT, lambda)
v_amb = lambda/( 4 * PRT );
% dv = lambda/(N_sweep * 4 * PRT);
if mod(N_sweep, 2) == 0
%     vel_axis = linspace(-v_amb, v_amb-dv, N_sweep);
    vel_axis = linspace(-N_sweep/2, N_sweep/2-1, N_sweep)./N_sweep .* 2 .* v_amb;
else
    vel_axis = linspace(-v_amb, v_amb, N_sweep);
end
dv = vel_axis(2) - vel_axis(1);
    for i = 1:length(r)
        for l = 1:length(phi)
            F = squeeze(Rr((l - 1)*N_sweep+1:l*N_sweep, i, l));
            F_fft = 1./sqrt(N_sweep) .* fftshift(fft(F));
%             figure; plot(vel_axis, abs(F_fft).^2);
            PT(i, l) = sum(abs(F_fft).^2 .* dv);
            Mu(i, l) = 1./PT(i, l) .* sum(vel_axis.' .* abs(F_fft).^2 .* dv);
            Sigma(i, l) = sqrt(1./PT(i, l) .* sum((vel_axis.' - Mu(i, l)).^2 .* abs(F_fft).^2 .* dv));
        end
    end
end