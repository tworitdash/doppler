function [a, r, itern] = MHradarOne(y, iter, data, t_avail, x0, sigma_n)
    obs = y;
    x = y(1);
    a = 0;
    r = 0;
    itern = 0;
    for i = 1:iter
        x_new = TMOne(x, obs(2));
        x_lik = LLOne(x, data, t_avail, x0, sigma_n);
        x_new_lik = LLOne(x_new, data, t_avail, x0, sigma_n);
        if (acceptance(real(x_lik) + log(priorOne(x, obs(1))), real(x_new_lik)+log(priorOne(x_new, obs(1)))))
            x = x_new;
            a = [a x_new];
            itern = [itern i];
        else
            r = [r x_new];
        end
    end
end