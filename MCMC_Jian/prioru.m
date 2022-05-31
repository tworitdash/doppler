function res = prioru(x, H, L)
    cond(1) = (x(1) > L(1)) && (x(1) < H(1));
    for i = 2:length(L)
        cond(i) = cond(i - 1) && (x(i) > L(i)) && (x(i) < H(i));
    end

    if cond(end)
        res = 1;
    else
        res = 0;
    end

end