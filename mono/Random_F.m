clear;

dx = 10;
dy = 10;

x = -5e3:dx:5e3;
y = -5e3:dy:5e3;
Lx = 1e3;
Ly = 0.5e3;
Corrx = exp(-abs(x)/Lx);
Corry = exp(-abs(y)/Ly);

figure; plot(x, Corrx); hold on; plot(y, Corry, 'o');