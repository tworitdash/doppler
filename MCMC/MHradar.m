function [accepted_umean, rejected_umean, accepted_usigma, rejected_usigma] = MHradar(x, iter, data, t_avail, x0, y0, N, sigma_n)
    obs = x;
    accepted_umean = 0;
    rejected_umean = 0;
    accepted_usigma = 0;
    rejected_usigma = 0;
    for i = 1:iter
        x_new = TMRadar(x);
        x_lik = LLRadar(x, data, t_avail, x0, y0, N, sigma_n);
        x_new_lik = LLRadar(x_new, data, t_avail, x0, y0, N, sigma_n);
        if (acceptance(real(x_lik) + log(priorRadar(x, obs)), real(x_new_lik)+log(priorRadar(x_new, obs))))
            x = x_new;
            accepted_umean = [accepted_umean x_new(1)];
            accepted_usigma = [accepted_usigma x_new(2)];
        else
            rejected_umean = [rejected_umean x_new(1)];
            rejected_usigma = [rejected_usigma x_new(2)];
        end
    end
end