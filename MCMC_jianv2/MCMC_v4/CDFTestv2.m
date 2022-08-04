clear;
N = 10000;
x = rand(1, N);
Mu = 0;
Sigma = 1;

U = normrnd(Mu, Sigma, [1 N]);

% F = 1/2 * (1 + erf((x - Mu)./(sqrt(2) .* Sigma)));

F = erfinv(2*x-1) * Sigma * sqrt(2) + Mu;
figure; histogram(U); hold on; histogram(F);
