function res = priorRadar(x, obs)

    if (x(2) > 0) && (x(2) < obs(2)) && (x(1) > (obs(1) - obs(2))) && (x(1) < (obs(1) + obs(2))) && (x(1) > 0)
        res = 1;
    else
        res = 0;
    end

end