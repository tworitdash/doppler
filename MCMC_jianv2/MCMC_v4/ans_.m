a = 1; % for example
n = 10:10:50;
x = a + (-1:.001:1).';
x(x==a) = a+eps; % avoid divide by zero
Z1 = csc((x-a)/2).*sin((n+1/2).*(x-a))/2/pi;
figure
plot(x,Z1)