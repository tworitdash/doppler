function [y, pynew, px] = TMu(x, Esigma)
%     y = [normrnd(x(1), obs(2), 1) normrnd(x(2), 5, 1) normrnd(x(3), 10, 1)];
    for i = 1:length(Esigma)
        y(i) = normrnd(x(i), Esigma(i), 1); % Draw a new value
        
        
        px(i) = normpdf(x(i), y(i), Esigma(i)); % Find the prior probability of x
        pynew(i) = normpdf(y(i), x(i), Esigma(i)); % Find the prior probability of y
    end 
end