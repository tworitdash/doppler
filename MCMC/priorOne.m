function res = priorOne(x, obs)
     if (x > obs-2) && (x < obs + 2)
        res = 1;
    else
        res = 0;
    end

end