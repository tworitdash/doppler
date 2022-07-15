K = 0:1:10;
n = 0:1:30;

for k = 1:length(K)
    
    F(k) = sum((-1).^n .* K(k).^(2.*n)./(factorial(n))); % + k.^4./2 - k.^6./6 + k.^8./24 - k.^10./120 + k.^12./720 - k.^14./5040;

end
figure; plot(K, db(F), '*'); hold on; plot(K, db(exp(-K.^2)));