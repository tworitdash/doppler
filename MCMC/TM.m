function y = TM(x, Esigma)
%     y = [normrnd(x(1), obs(2), 1) normrnd(x(2), 5, 1) normrnd(x(3), 10, 1)];
    for i = 1:length(Esigma)
        y(i) = normrnd(x(i), Esigma(i), 1);
    end 
end