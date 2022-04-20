function [accepted, rejected, itern, E] = MHu(E, iter, data, t_avail, x0, sigma_n)

%% Inputs:

% E - Struct of MCMC options 
% iter - number of iterations
% data - data available
% t_avail - time instances of the data
% x0 - start position of the target - assumed to be known
% sigma_n - noise standard deviation - assumed to be known

%% Outputs:

% accepted - accpetd values
% rejected - rejected values 
% itern - iterations at which we accept a value
% E - the new Struct of MCMC options after processing - A new E.sig is
% updated


    x = E.E0; % save the initial value in x
    accepted = zeros(1, E.n); 
    rejected = zeros(1, E.n);
    
    
    itern = 0;
    
%     fig1 = figure('Position',[0 0 800 800]);
%     fig2 = figure('Position',[850 0 800 800]);
%     
%     figure(fig1); histogram(normrnd(E.gt(1), E.gt(2), [1 5000])); 
%     figure(fig2); histogram(normrnd(E.gt(3), E.gt(4), [1 5000]));
    
    for i = 1:iter
        disp(i);
        [x_new, pxnew, px] = TMu(x, E.sig);     % x_new is a new sample drawn for u 
                                                % pxnew is the prior
                                                % prior probability of x_new
                                                % prior probability of x 
        
        x_lik = LLu(x, data, t_avail, x0, sigma_n, px); % Likelihood of x wth data
        x_new_lik = LLu(x_new, data, t_avail, x0, sigma_n, pxnew); % Likelihood of x_new with data
        
%         if mod(i, 10) == 0
%             E.sig = 0.95 .* E.sig;
%         end
        
%         if (acceptance(real(x_lik) + log(prioru((x), E.H, E.L)), real(x_new_lik)+log(prioru((x_new), E.H, E.L))))...
%                 && (acceptance(imag(x_lik) + log(prioru((x), E.H, E.L)), imag(x_new_lik)+log(prioru((x_new), E.H, E.L))))
            
   
%% Acceptance Logic

%         if (acceptance(real(x_lik), real(x_new_lik))) || (acceptance(imag(x_lik), imag(x_new_lik)))
        if (acceptance((x_lik), (x_new_lik)))
            
            accepted = [accepted; x_new];
            itern = [itern i];
%             
%             if length(itern) > 2
%                 E.sig = abs(x_new - x);
%             end
            x = x_new;
            
           
            
%             if length(itern) > 1
%                 ep = abs(itern(end) - itern(end -1));
%                 if ep < 1
%                     disp('conv');
%                     E.sig = 0.2;
%                 end
%             end
            
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