function [y, pynew, px, Esig] = TMu(x, Esigma)
%     y = [normrnd(x(1), obs(2), 1) normrnd(x(2), 5, 1) normrnd(x(3), 10, 1)];
    for i = 1:length(Esigma)
        
        if i < 3
            y(i) = Esigma(i) .* rand();
            px(i) = 1/Esigma(i);
            pynew(i) = 1/Esigma(i);
            
        else
        
            Esigma(i) = 200 * betarnd(y(1), y(2), 1);
            Esig = Esigma(i);

            y(i) = normrnd(x(i), Esigma(i), 1); % Draw a new value
            px(i) = normpdf(x(i), x(i), Esigma(i)); % Find the prior probability of x
            pynew(i) = normpdf(y(i), x(i), Esigma(i)); % Find the prior probability of y
        end
    end
end