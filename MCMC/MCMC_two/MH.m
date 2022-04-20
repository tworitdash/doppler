function [accepted, rejected, itern] = MH(E, iter, data, t_avail, x0)
    x = E.E0;
    accepted = zeros(1, E.n);
    rejected = zeros(1, E.n);
    itern = 0;
    
    fig1 = figure('Position',[0 0 800 800]);
    fig2 = figure('Position',[850 0 800 800]);
    
%     figure(fig1); histogram(normrnd(E.gt(1), E.gt(2), [1 5000])); 
%     figure(fig2); histogram(normrnd(E.gt(3), E.gt(4), [1 5000]));
    
    for i = 1:iter
        x_new = TM(x, E.sig);
        x_lik = LLsig(x, data, t_avail, x0);
        x_new_lik = LLsig(x_new, data, t_avail, x0);
        if (acceptance(real(x_lik) + log(prior(x(1:2), E.L(1:2), E.H(1:2))), real(x_new_lik)+log(prior(real(x_new), E.L(1:2), E.H(1:2)))))...
                && (acceptance(imag(x_lik) + log(prior(x(3:4), E.L(3:4), E.H(3:4))), imag(x_new_lik)+log(prior(x_new(3:4), E.L(3:4), E.H(3:4)))))
            x = x_new;
            accepted = [accepted; x_new];
            itern = [itern i];
            
            figure(fig1);
            
%             hold on; 
            histogram(normrnd(x_new(1), x_new(2), [1 5000])); 
            hold off;
            
            figure(fig2); 
%             hold on; 
            histogram(normrnd(x_new(3), x_new(4), [1 5000])); 
            hold off;
            
        else
            rejected = [rejected; x_new];
        end
    end
end