function [unew, punew, pu] = TMu(u, Esigma)

    for i = 1:length(Esigma)
        unew(i) = normrnd(u(i), 0.0001, 1); % Draw a new value
        
        
%         pu(i) = normpdf(u(i), u(i), Esigma(i)); % Find the prior probability of x
%         punew(i) = normpdf(unew(i), u(i), Esigma(i)); % Find the prior probability of y
        
%         pu(i) = normpdf(u(i), 0, Esigma(i)); % Find the prior probability of x
%         punew(i) = normpdf(unew(i), 0, Esigma(i)); % Find the prior probability of y
        pu(i) = -0.5*log(2*pi)-0.5*(u(i))^2;
        punew(i) = -0.5*log(2*pi)-0.5*(unew(i))^2;
        
    end 
end