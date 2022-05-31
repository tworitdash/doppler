function [y, pynew, px] = TMu_uniformprior(E)
%     y = [normrnd(x(1), obs(2), 1) normrnd(x(2), 5, 1) normrnd(x(3), 10, 1)];
    for i = 1:length(E.L)
        y(i) = E.L + (E.H - E.L) .* rand(); % Draw a new value
        
        
%         px(i) = normpdf(x(i), y(i), Esigma(i)); % Find the prior probability of x
%         pynew(i) = normpdf(y(i), x(i), Esigma(i)); % Find the prior probability of y
        px(i) = 1; pynew(i) = 1;
    end 
end