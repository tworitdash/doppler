function ret = acceptance_0(x, x_new)
   
        alpha = exp(x_new - x);
        rnd = rand;
        
    if alpha > rnd
        ret = 1;
    else
        ret = 0;
        %accept = rand;
        %ret = (accept < exp(x_new - x));
%         ret = 0;
    end
end