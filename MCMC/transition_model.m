function y = transition_model(x)
    y = [x(1) normrnd(x(2), 0.5, 1)];
end