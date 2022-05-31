function [accepted, rejected, itern, E, sample_MC_seq, sample_all, Flag] = MHu_uniformpriordouble(E, iter, data, t_avail, x0, sigma_n, epislon)

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
    sample_MC_seq = zeros(E.n, iter);
    sample_all = zeros(E.n, iter);
    
    itern = [];
    f = 1;
    Flag = 0;
    
%     fig1 = figure('Position',[0 0 800 800]);
%     fig2 = figure('Position',[850 0 800 800]);
%     
%     figure(fig1); histogram(normrnd(E.gt(1), E.gt(2), [1 5000])); 
%     figure(fig2); histogram(normrnd(E.gt(3), E.gt(4), [1 5000]));
    
    for i = 1:iter
        if Flag == 0
            f = 1;
        end
        
        if mod(i,5000)==0
        disp(i);
        end
        [x_new, pxnew, px] = TMu_uniformpriordouble(E);     % x_new is a new sample drawn for u 
                                                % pxnew is the prior probability of x_new
                                                
                                                % px is the prior probability of x 
         for l = 1:E.n 
            sample_all(l, i) = x_new(l);
         end
         
%         sample_all(i) = x_new;

        x_lik = LLudouble(x, data, t_avail, x0, sigma_n, pxnew); % Likelihood of x wth data
        x_new_lik = LLudouble(x_new, data, t_avail, x0, sigma_n, px); % Likelihood of x_new with data
        
%         if mod(i, 10) == 0
%             E.sig = 0.95 .* E.sig;
%         end
        
%         if (acceptance(real(x_lik) + log(prioru((x), E.H, E.L)), real(x_new_lik)+log(prioru((x_new), E.H, E.L))))...
%                 && (acceptance(imag(x_lik) + log(prioru((x), E.H, E.L)), imag(x_new_lik)+log(prioru((x_new), E.H, E.L))))
            
   
%% Acceptance Logic

%         if (acceptance(real(x_lik), real(x_new_lik))) || (acceptance(imag(x_lik), imag(x_new_lik)))
        if (acceptance(x_lik, x_new_lik)) itern = [itern i];
            
            for l = 1:E.n
%                 accepted(l, :) = [accepted(l, :); x_new(l)];
                accepted(l, i) = [x_new(l)];
                sample_MC_seq(l, i) = x_new(l);
                
            end

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
        for l = 1:E.n
%             rejected(l, :) = [rejected(l, :); x_new(l)];
            rejected(l, i) = [x_new(l)];
            sample_MC_seq(l, i) = x(l);
        end
        end
        
%         if i > 1
%             slope(i-1) = sample_MC_seq(i) - sample_MC_seq(i -1);
%             if i > 2100
%                 while f
%                     if slope(i-2000:i-1) < epislon
%                         Flag = 1;
%                         NB = round(mean([i i - 2000]));
%                         f = 0;
%                         
%                     else
%                         NB = 0;
%                         Flag = 0;
%                         f = 0;
%                     end
%                 end
%             end
%         end
    end
end