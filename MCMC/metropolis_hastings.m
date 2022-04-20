function [accepted, rejected] = metropolis_hastings(x, iter, data)
    accepted = 0;
    rejected = 0;
    for i = 1:iter
        x_new = transition_model(x);
        x_lik = manual_log_like_normal(x, data);
        x_new_lik = manual_log_like_normal(x_new, data);
        if (acceptance(x_lik + log(prior(x)), x_new_lik+log(prior(x_new))))
            x = x_new;
            accepted = [accepted x_new(2)];
        else
            rejected = [rejected x_new(2)];
        end
    end
end