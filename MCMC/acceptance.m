function ret = acceptance(x, x_new)
    if x_new > x
        ret = 1;
    else
        accept = rand;
        ret = (accept < exp(x_new - x));
    end
end