clear;


N = linspace(1, 257, 257);


for i = 1:length(N)
    
    x = linspace(1, N(i), N(i));

    calc_fft = @() fft(x, N(i));

    time_consumed(i) = timeit(calc_fft);
end

plot(N, time_consumed); grid on;