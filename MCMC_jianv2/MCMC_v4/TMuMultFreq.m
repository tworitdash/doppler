function [y, pynew, px] = TMuMultFreq(E, theta)
%     y = [normrnd(x(1), obs(2), 1) normrnd(x(2), 5, 1) normrnd(x(3), 10, 1)];
        pd = makedist('Poisson', 'lambda', theta(1));
        t = truncate(pd, 1, 50);
        y(1) = random(t, 1);
            
        px(1) = poisspdf(theta(1), theta(1));
        pynew(1) = poisspdf(y(1), theta(1));
        
    for i = 2:y(1)+1
    
            y(i) = E.L(2) + (E.H(2) - E.L(2)) .* rand; % Draw a new value


    %         px(i) = normpdf(x(i), y(i), Esigma(i)); % Find the prior probability of x
            px(i) = 1;
            pynew(i) = 1;
    %         pynew(i) = normpdf(y(i), x(i), Esigma(i)); % Find the prior probability of y
       
     end 
end