function [accepted, rejected, itern] = MHu(E, iter, data, t_avail, x0)
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
        if (acceptance(real(x_lik) + log(prior(real(x, E.H, E.L))), real(x_new_lik)+log(prior(real(x_new, E.H, E.L)))))...
                && (acceptance(imag(x_lik) + log(prior(imag(x, E.H, E.L))), imag(x_new_lik)+log(prior(imag(x_new, E.H, E.L)))))
            x = x_new;
            accepted = [accepted; x_new];
            itern = [itern i];
            
%             figure(fig1);
%             
% %             hold on; 
%             histogram(normrnd(x_new(1), x_new(2), [1 5000])); 
%             hold off;
            
%             figure(fig2); 
% %             hold on; 
%             histogram(normrnd(x_new(3), x_new(4), [1 5000])); 
%             hold off;
%             
        else
            rejected = [rejected; x_new];
        end
    end
end