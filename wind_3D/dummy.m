T = 1;

fs = 1024;

t = linspace(0, T-1/fs, fs);



% y = zeros(size(t));
% 
% y(30:60) = 1;

c1=[1 -1 1 -1];

Npoint=size(t,2)/4;

c2 = ones(Npoint,1).';

c3 = kron(c1,c2);

c4= (c3+1).*(pi/2);

y=c4;

plot(t, y); grid on;

N = length(t);

fmax = N .* fs;

yf = 1/sqrt(N) .* fft(fftshift(y));

f = linspace(-fmax/2, fmax/2-fs, length(t));

k = 4e10;

yf_new = yf .* exp(-1j .* pi .* f.^2 ./k);

y_new = ifft(fftshift(sqrt(N) .* yf_new));

hold on;
plot(t, abs(y_new));

%% 

a =  exp(-1j .* pi .* f.^2 ./k);

at = ifft(fftshift(sqrt(N) .* a));


bt = exp(-(k .* t.^2)./(4 .* 1j .* pi))./(sqrt(1j/k) .* sqrt(2.*pi));
% 
% btn = conv(bt, y, 'same');

figure; plot(t, abs(at)); % hold on; %plot(t, abs(bt)); hold on; plot(t, abs(btn));
figure; plot(t, abs(bt))

