clear;
close all;


m_ = load('markers.mat'); 
markers = m_.markers;
% N = 2.^(linspace(2, 6));

N = 5:65;
sigma_true = linspace(eps, 7.5, 100);
v_amb = 7.5;

for i = 1:length(N)
    for k = 1:length(sigma_true)

        PRT = 1e-3;
        lambda = 0.03;
        t = 0:PRT:N(i)*PRT;
        MC = 128;
        
        for m = 1:MC

            mu = normrnd(0, sigma_true(k), [1 length(t)]);
            n_sig = 0; 
    %         n_sig = sqrt(0.000001);

            [s] =  TD_generatorV2(mu, lambda, t, n_sig);

            if (mod(N(i), 2) == 0)
                    hsnew = N(i) +1;
                    %     vel_axis_hs = linspace(-hs/2, hs/2-1, hs)./hs .* 2 .* v_amb;
                    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
            else
                    hsnew = N(i);
                    vel_axis_hs = linspace(-v_amb, v_amb, hsnew);
            end

            S = 1./sqrt(hsnew) .* fftshift(fft(squeeze(s), hsnew));


            diff_v = diff(vel_axis_hs); dv = diff_v(1); dv_(i) = dv;
            sigmabydv(k, i) = sigma_true(k)/dv_(i);

            PT_i = squeeze(abs(S)).^2;
            PT(m, i, k) = sum(PT_i .* dv);
            mu_re_i = vel_axis_hs .* squeeze(abs(S)).^2;
            mu_re_(m, i, k) = sum(mu_re_i .* dv)./PT(m, i, k);
            sigma_re_i = (vel_axis_hs - mu_re_(m, i, k)).^2 .* squeeze(abs(S)).^2 .* dv;
            sigma_re_(m, i, k) = sqrt(sum(sigma_re_i)./PT(m, i, k));
        end
        mu_re(k, i) = squeeze(mean(mu_re_(:, i, k), 1));
        sigma_re(k, i) = squeeze(mean(sigma_re_(:, i, k), 1));
        
        sigma_error(k, i) = sigma_re(k, i)./sigma_true(k);
    end
end
%%
% txt = ['Retrieved \sigma'];
% 
% xl = 'Number of points N';
% yl = 'dv/\sigma';
% zl = '\sigma_{retrieved}';
% % sigma_error = abs(sigma_re.'- (repmat(sigma_true, length(N), 1)).')./((repmat(sigma_true, length(N), 1)).') * 100;
% 
% N_ = repmat(N, length(sigma_true), 1);
% sigma_true_ = repmat(sigma_true, length(N), 1);
% 
% 
% % surplot_pcolor(N_, sigma_true_.', sigma_error, xl, yl, zl, txt);
% surplot_pcolor(N_, sigmabydv, sigma_error, xl, yl, zl, txt); 
% caxis([min(min(sigma_error)) max(max(sigma_error))]);
% ylim([0 8]);

%% 1D plot with respect to the sigma/dv at a fixed N

n_idx = 5;

sigma_error_1D = sigma_error(:, n_idx);
sigmabydv_1D = sigmabydv(:, n_idx);

txt = ['Percentage error in \sigma estimation'];

dtext = ['N = ', num2str(N(n_idx))];
xl = '\sigma/dv';

f = figure(106); hold on; f.Position = [10 10 1000 1000];
color = 'k';
yl =  ['\sigma_{error} in % at \phi = \phi_{wind}'];
% plott(sigmabydv_1D, sigma_error_1D, xl, yl, txt, 2, dtext, color, markers(2))
plott(sigmabydv_1D, sigma_re(:, n_idx), xl, yl, txt, 2, dtext, color, markers(3))

% xlim([0 8]);
