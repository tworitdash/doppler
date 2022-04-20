function [am, rm, as, rs, an, rn itern] = MHradarOne(y, iter, data, t_avail, x0, y0, sigma_n)
    obs = y;
    x = y;
    am = 0;
    rm = 0;
    as = 0;
    rs = 0;
    an = 0; 
    rn = 0;
    itern = 0;
    for i = 1:iter
        x_new = TMOne(x, obs);
        x_lik = LLOne(x, data, t_avail, x0, y0, sigma_n);
        x_new_lik = LLOne(x_new, data, t_avail, x0, y0, sigma_n);
        if (acceptance(real(x_lik) + log(priorOne(x, obs)), real(x_new_lik)+log(priorOne(x_new, obs))))
            x = x_new;
            am = [am x_new(1)];
            as = [as x_new(2)];
            an = [an x_new(3)];
            itern = [itern i];
        else
            rm = [rm x_new(1)];
            rs = [rs x_new(2)];
            rn = [rn x_new(3)];
        end
    end
end