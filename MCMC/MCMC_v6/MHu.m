function [accepted, rejected, itern, E] = MHu(E, iter, data, t_avail, r0, sigma_n)

%% Inputs:

% E - Struct of MCMC options 
% iter - number of iterations
% data - data available
% t_avail - t instances of the data
% r0 - start position of the target - assumed to be known
% sigma_n - noise standard deviation - assumed to be known

%% Outputs:

% accepted - accpetd values
% rejected - rejected values 
% itern - iterations at which we accept a value
% E - the new Struct of MCMC options after processing - A new E.sig is
% updated


    u = E.E0; % save the initial value in x
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
        [u_new, punew, pu] = TMu(u, E.sig);     % x_new is a new sample drawn for u 
                                                % prior probability of x_new
                                                % prior probability of x 
        
        u_lik = LLu(u, data, t_avail, r0, sigma_n, pu); % Likelihood of x wth data
        u_new_lik = LLu(u_new, data, t_avail, r0, sigma_n, punew); % Likelihood of x_new with data
   
%% Acceptance Logic

        if (acceptance((u_lik), (u_new_lik)))
            
            
            accepted = [accepted; u_new];
            itern = [itern i];             
            u = u_new;
        else
            rejected = [rejected; u_new];
        end
    end
end