clear;

lambda = 3e-2;


T = 3e-3;
t = [eps:1e-3:T] ; % 0.1:1e-3:0.1+T];

u = 3;

a = 4*pi/lambda * u;

x = a + 4*pi/lambda * linspace(-15, 15, 1000);

x(x==a) = a+eps;

for i = 1:length(x)

    Z(i) = sum(-((cos(a.*t) - cos(x(i).*t)).^2 + (sin(a.*t) - sin(x(i).*t)).^2));

end

% Z1 = csc((a+x)/2) .* sin(1/2.*(2.*T + 1).*(a + x)) -1/2 .*  csc(a) .* sin(2.*a.*T+a) - 1/2 .* csc(x) .* sin(2.*T.*x+x)
Z1 = csc((x-a)/2) .* sin((T+1/2).*(x-a)); %/2/pi;%- 2*T -1;

figure; plot(x*lambda/(4 * pi), (Z)); 
hold on; plot(x*lambda/(4 * pi), (Z1));

% n = 100;
% 
% Z = 1/(2*pi) * csc((x - u)/2) .* sin((n+1/2).*(x-u));
% 
% figure; plot(x, Z);