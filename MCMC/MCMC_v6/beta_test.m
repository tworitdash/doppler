%% Beta distribution test

for i = 1:100000

alpha = rand() * 5;
beta = rand() * 5;

sig(i) = 7 * betarnd(alpha, beta, 1);

end

figure; plot(sig); figure; histogram(sig);