function [accepted, rejected, itern, E, sample_MC_seq, sample_all] = MHMultFreq(E, iter, data, t_avail, x0, y0, sigma_n)

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


    theta = E.E0; % save the initial value in x
    accepted = [];%zeros(1, E.n); 
    rejected = [];%zeros(1, E.n);
%     sample_MC_seq = zeros(E.n, iter);
%     sample_all = zeros(E.n, iter);
    
    itern = [];
    
    for i = 1:iter
        
        if mod(i,5000)==0
            disp(i);
        end
        
        [theta_new, pthetanew, ptheta] = TMuMultFreq(E, theta);     % x_new is a new sample drawn for u 
                                                % pxnew is the prior probability of x_new
                                                
                                                % px is the prior probability of x 

        sample_all(i).seq = theta_new;

        theta_lik = LLMultFreq(theta, data, t_avail, x0, y0, sigma_n, pthetanew); % Likelihood of x wth data
        theta_new_lik = LLMultFreq(theta_new, data, t_avail, x0, y0, sigma_n, ptheta); % Likelihood of x_new with data
       
   
%% Acceptance Logic

%         if (acceptance(real(x_lik), real(x_new_lik))) || (acceptance(imag(x_lik), imag(x_new_lik)))
        if (acceptance(theta_lik, theta_new_lik))         
            
%             accepted = [accepted; theta_new];
            itern = [itern i];
            sample_MC_seq(i).seq = theta_new;
            theta = theta_new;
            
%             
        else
%             rejected = [rejected; theta_new];
            sample_MC_seq(i).seq = theta;
        end
    end
end