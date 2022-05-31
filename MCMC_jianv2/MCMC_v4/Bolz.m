%% Testing the Gibbs/ Bolzmann Distribution

y = [1 2 3 2 1 3 4 5 7 8];
N = length(y);


data = y + 1 * randn(1, N);


% figure; plot(y); hold on; plot(data);

c = 10;
C = linspace(1, c, c);

beta = linspace(1, 1000, 100);

for k = 1:10

for i = 1:1
    
yj = normrnd(y(1), 2, c);

Ej = abs(data(i) - y).^2;


PrxC = (exp(-beta(k) * Ej)+eps)/(sum(exp(-beta(k)*Ej))+eps);

hold on; plot(C, (PrxC)); grid on;

end
end