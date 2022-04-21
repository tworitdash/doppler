function [unew, punew, pu] = TMu(u, Esigma)

    for i = 1:length(Esigma)
        unew(i) = normrnd(u(i), Esigma(i), 1); % Draw a new value
        
        
        pu(i) = normpdf(u(i), u(i), Esigma(i)); % Find the prior probability of x
        punew(i) = normpdf(unew(i), u(i), Esigma(i)); % Find the prior probability of y
    end 
end