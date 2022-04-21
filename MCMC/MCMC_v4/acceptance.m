function ret = acceptance(u, u_new)
   
        alpha = exp(u_new - u);
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