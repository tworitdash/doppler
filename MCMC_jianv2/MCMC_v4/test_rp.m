clear; 
close all;



Nu = 100000;

% histogram(x);

MC = 1;

for k = 1:100000

for i = 1:MC
    x = -pi + 2 * pi * rand(1, k);
    y = exp(1j .* x); 
    s(i) = sum(y);
end
    Sm(k) = mean(s);
    Sstd(k) = std(s);
end

% y =  sin(x);

%/Nu;

%% 

% z = normrnd(5, 0.1, [1 Nu]);
% 
% c = sum(z)/Nu;

