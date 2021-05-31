from DS_simulatorV2 import DS_simulatorV2R

SNR = 10
m0 = 1
mu = 5
sigma = 0.2
n = 1024
v_amb = 7.5
N = 128

[data, data_f, data_f_full, data_full, X_full, Theta_full] = DS_simulatorV2(SNR, m0, mu, sigma, n, v_amb, N)