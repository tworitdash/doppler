function [mu, sigma, mu_r_sigma, sigma_r_sigma, r_new, dv] = retrieval_V2(s, v_amb, hs, Nr, Nphi, sec, r_avg_op, r, Ravg, multi_rot_op, mu_prior, phi_axis)

if (mod(hs, 2) == 0)
    hsnew = hs+1;
%     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
else
    hsnew = hs;
    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
end
n_rot = size(s, 1);

S1_norm = zeros(n_rot, Nr, Nphi, hsnew);

for n = 1:n_rot
    for i = 1:Nr
        si = squeeze(s(n, i, :));
        si2 = reshape(si, [hs Nphi]).';
        z = zeros(Nphi, hsnew-hs);
        si3 = [si2 z];
        si4 = (reshape(si3.', [1 Nphi*hsnew]));

        [S1, F1, Ti1, P1] = spectrogram(si4, hsnew, 0, hsnew, hsnew*sec*1/60); 
        S1_norm(n, i, :, :) = 1./sqrt(hsnew) .* fftshift(S1',2);
%         S1_norm(n, i, :, :) = S1f(i, :, :)./max(max(squeeze(S1f(i, :, :))));
        S1_norm_db = 20*log10(abs(squeeze(S1_norm(n, i, :, :))));

    %     if (n == 1) || (n == length(Vmean))
    % 
        txt = ['Spectrum']; %, 'Vmean = ', num2str(Vmean(n))];
    
        xl = 'Doppler Velocity [m.s^{-1}]';
        yl = 'Azimuthal Angle \phi [{\circ}]';
        zl = 'Power [dB]';
    
        surplot_pcolor(vel_axis_hs, phi_axis*180/pi, S1_norm_db, xl, yl, zl, txt);
    %     
    %     end
    end
end


diff_v = diff(vel_axis_hs); dv = diff_v(1);
% PT = zeros(Nr, length(phi_axis));
% mu_re = zeros(Nr, length(phi_axis));
% sigma_re = zeros(1, length(phi_axis));
for n = 1:n_rot
    for i = 1:Nr
        for k = 1:Nphi
            vel_axis_hs_ = mean(mu_prior(n, i, (k - 1)*hs+1:k*hs))+vel_axis_hs;
%             vel_axis_hs_ = vel_axis_hs;
            PT_i = squeeze(abs(S1_norm(n, i, k, :))).^2;
            PT(n, i, k) = sum(PT_i .* dv);
            mu_re_i = vel_axis_hs_.' .* squeeze(abs(S1_norm(n, i, k, :))).^2;
            mu_re(n, i, k) = sum(mu_re_i .* dv)./PT(n, i, k);
            sigma_re_i = (vel_axis_hs_.' - mu_re(n, i, k)).^2 .* squeeze(abs(S1_norm(n, i, k, :))).^2 .* dv;
            sigma_re(n, i, k) = sqrt(sum(sigma_re_i)./PT(n, i, k));
        end
    end
end

if multi_rot_op == 1 
    Mu_re_avg = squeeze(mean(mu_re, 1));
    Sigma_re_avg = squeeze(mean(sigma_re, 1));
    Mu_re_std_avg = std(mu_re, 1, 1);
    Sigma_re_std_avg = std(sigma_re, 1, 1);
else
    Mu_re_avg = squeeze(mu_re(1, :, :));
    Sigma_re_avg = squeeze(sigma_re(1, :, :));
    Mu_re_std_avg = zeros(Nr, Nphi);
    Sigma_re_std_avg = zeros(Nr, Nphi);
end


if r_avg_op == 1
    [~, numr] = find(r < Ravg);
    Num_N = length(numr);
    r_new = r(1:Num_N:end);
    
    for rng = 1:length(r_new)-1
        Mu_re_rng(rng, :) = mean(Mu_re_avg((rng - 1)*Num_N+1:rng*Num_N, :),  1);
        Sigma_re_rng(rng, :) = mean(Sigma_re_avg((rng - 1)*Num_N+1:rng*Num_N, :),  1);
        Mu_re_std_rng(rng, :) = std(Mu_re_avg((rng - 1)*Num_N+1:rng*Num_N, :), 1, 1);
        Sigma_re_std_rng(rng, :) = std(Sigma_re_avg((rng - 1)*Num_N+1:rng*Num_N, :), 1, 1);
    end
    
    mu = Mu_re_rng;
    sigma = Sigma_re_rng;
    mu_r_sigma =  Mu_re_std_rng;
    sigma_r_sigma = Sigma_re_std_rng;
    
    
else
    mu = Mu_re_avg;
    sigma = Sigma_re_avg;
    mu_r_sigma =  Mu_re_std_avg;
    sigma_r_sigma = Sigma_re_std_avg;
    r_new = r;
end
