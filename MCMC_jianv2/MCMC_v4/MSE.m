clear;
close all;


Mu = 7.5 .* 0.4189;
Sigma = 2 .*  0.4189;

Ntvec = 128;

% rng(1, 'twister');
Nu = 1000000;

for n = 1:length(Ntvec)
    
Nt = Ntvec(n);
K = 0:1:Nt-1;


x = zeros(1, Nt);

% U = normrnd(Mu, Sigma, [1 1]);

for m = 1:Nt
    beta = 2 * pi .* rand([1 Nu]); 
%     U = Mu - Sigma/2 + Sigma .* rand([1 Nu]);
    U = normrnd(Mu, Sigma, [1 Nu]);
    if m == 2
     figure(100); hold on; histogram(U + beta/K(m));% hold on; histogram(U);
    end
    x(m) = [sum(exp(1j .* K(m) .* U + 1j .* beta))];

end

% figure(1); hold on; plot(real(x)); hold on; plot(imag(x)); grid on;
% 
figure(2); hold on; histogram(real(x), 10);
figure(3); hold on; histogram(imag(x), 10);

Mr(n) = mean(real(x));
Sr(n) = std(real(x));

Mi(n) = mean(imag(x));
Si(n) = std(imag(x));



end

% figure; plot(Ntvec, Mr);
% figure; plot(Ntvec, Sr);

% figure(2); hold on; plot(db(abs(fft(x))));

%% 

X = linspace(-1, 24, 10000);
k_ = 1;

Y = Nu .* 1/(4*pi) * (erf((Mu - X + 2 * pi/K(2))./(sqrt(2) * Sigma)) ...
    - erf((Mu - X)./(sqrt(2) * Sigma)));

figure(100); hold on; plot(X, Y);

% sum(cos(X) .* Y)
% sum(sin(X) .* Y)


%% 

u = linspace(0, 2 * pi, 10000);
for m = 1:Nt
    if mod(K(m), 10) == 0
     figure(1000); hold on; plot(u, 1/sqrt(2*pi*Sigma^2) .* exp(-(K(m).*u - K(m).* Mu).^2/(2*K(m).^2 * Sigma.^2)));
    end
end
