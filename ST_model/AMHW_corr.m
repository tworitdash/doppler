function [Corr] = AMHW_corr(Theta, delta, tau)
    bs      = Theta(1);
    cs      = Theta(2);
    bt      = Theta(3);
    ct      = Theta(4);
    theta   = Theta(5);
    
    Corr_num    = exp(-(delta/bs).^cs) .* exp(-(tau/bt).^ct);
    Corr_denum  = 1 - theta .* (1 - exp(-(delta/bs).^cs)) .* (1 - exp(-(tau/bt).^ct));
    
    Corr        = Corr_num./Corr_denum;
end