N = 128;
[data, data_f, X, Theta] = DS_simulatorV3(10^(30/10), 1, 0, 0.2, N, 7.5);

df = 1/sqrt(N) * fftshift(fft(data));
df = df./max(abs(df));

figure; plot(linspace(-7.5, 7.5, N), db(abs(df))); grid on;