%% HD signal generator

close all;
clear;

SNR_db = 20;                % Noise is added with a given SNR in dB
SNR = 10^(SNR_db/10);       % SNR in linear scale


lambda = 1;             

r0 = 0;                     
u = 0.01;                      % Ground truth u

N = 100; 
M = 5;

K = linspace(1, 10, 10);

Nt = K(end) * N;                 % Number of Truth samples (N * M * 20)

dT = 0.1;                  % t step
r(1) = r0;                  

Z(1) = exp(1j * 4 * pi/lambda .* r(1)); 

t_whole = linspace(0, dT*(Nt - 1), Nt);



for i = 2:Nt
%     r(i) = r(i - 1) + u * dT;
    r(i) = r(1) + u * t_whole(i);
    Z(i) = exp(1j * 4 * pi/lambda .* r(i)); % Ground truth samples
end

Noise = sum(abs(Z).^2)./(Nt .* SNR);         % Finding Noise power and ...
                                             %noise variance with the data and given SNR
sigma_n = sqrt(Noise);

Z_model = Z + sigma_n .* (randn(1, Nt) + 1j .* randn(1, Nt))./sqrt(2); % Adding complex noise 


%% Available samples [Measurement model with only a few samples]

  


Z_model_re = reshape(Z_model, [N K(end)]); 
Z_avail = Z_model_re(1:M, :); 

Z_avail_vec = reshape(Z_avail, [M * K(end) 1]); % available samples for measurements


% Z_avail_vec = Z_avail_vec_ + sigma_n .* (randn(1, length(Z_avail_vec_)).'+ 1j .* randn(1, length(Z_avail_vec_)).')./sqrt(2);

for k = 1:K(end)
    t(:, k) = (k - 1) * N + [0:M-1]; % This for loop calculates the t instances of the available samples
end

t_avail = reshape(t, [length(Z_avail_vec) 1]) .* dT; % vectorize the available time instances

%% 

u_test = linspace(-0.05, 0.2, 200);
% u_test = 0.01;


for i = 1:length(u_test)
%     LL(i) = LLu(u_test(i), Z_avail_vec, t_avail, r0, sigma_n, 1, Z, t_whole, t_avail);
    pu = normpdf(u, 0, 1);
    LL(i) = LLu(u_test(i), Z_avail_vec, t_avail, r0, sigma_n, pu);
end

figure; plot(u_test, LL)



