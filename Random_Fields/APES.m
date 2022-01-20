function [H, beta, omega, g_mtr] = APES(sig, M, N_omg)
% Amplitude and Phase Estimation Method (APES)
% Input Parameters:
% sig --- the input signal sequence
% M --- a parameter to compute averaged convariance matrix
 
% N_omg --- the number of frequency samples in the frequency domain
% 
% Output:
% H --- the computed filter for APES
 
% beta --- the vector of signal amplitudes of different frequency components
 
% omega --- the vector of angular frequencies
% g_mtr --- matrix formed by the vectors of g obtained with different
% signal segments ?intermediate result used for correctness check?
%
% Reference: P. Stoica , R. L. Moses, Spectral Analysis of Signals, Pearson
% Prentice Hall, 2005.
%
% Author: Jianping Wang
% Date: Jan 31, 2021
%
ii = sqrt(-1);
if size(sig,1) == 1
sig = sig.';
end
N = length(sig);
R = zeros(M,M);
sig_mtr = flipud(hankel(sig(1:M), sig(M:N)));
for k = 1:N-M+1
R = R + sig_mtr(:,k) * sig_mtr(:,k)';
end
R = R/(N-M+1);
%R_inv = inv(R);
omega = 2*pi/N_omg*linspace(0, N_omg-1, N_omg);
g_mtr = zeros(M, N_omg);
a = exp(-ii*2*pi/M*(0:M-1)).' ;
H = zeros(M,N_omg);
beta = zeros(N_omg,1);
for k = 1:N_omg
omega_ele = omega(k);
exp_omega = diag(exp(-ii * omega_ele * linspace(M,N, N-M+1) ));
g_vec = 1/(N-M+1) * sum(sig_mtr * exp_omega, 2);
g_mtr(:,k) = g_vec;
temp = 1 - g_vec' * (R \ g_vec);
term1 = temp * (R \ a);
term2 = (g_vec' * (R \ a)) * (R \ g_vec);
term3 = temp * a' * (R \ a);
term4 = a' * (R \ g_vec);
H(:,k) = (term1 + term2)/(term3 + abs(term4)^2 );
beta(k) = term4/(term3 + abs(term4)^2 );
end
end

