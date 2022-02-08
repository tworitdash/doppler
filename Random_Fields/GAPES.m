function [H, beta, sig_out] = GAPES(sig, I_a, M, N_omg)
% Amplitude and Phase Estimation Method with Gapped Data (GAPES)
% Input Parameters:
% sig --- the input signal sequence
% M --- a parameter to compute averaged convariance matrix
% I_a --- indices of available samples
% N_omg --- Number of angular frequencies
%
% Output:
% H --- the computed filter
% beta --- the vector of signal amplitudes of different frequency
% components
% sig_out --- full signal vector (composed of signal samples available and recovered samples)
%
% Reference: P. Stoica , R. L. Moses, Spectral Analysis of Signals, Pearson
% Prentice Hall, 2005.
%
% Author: Jianping Wang
% Date: Jan 31, 2021
%
ii = sqrt(-1);
if size(sig,1)==1
sig = sig.';
end
N_s = length(sig);
%====================================================
% Either of the two lines below can be used
[H, beta, omega, ~] = APES(sig, M, N_omg);
% [H, beta, omega, ~] = APES_FFT(sig, M, N_omg);
%====================================================
err_i0 = error_APES(sig, H, beta, omega);
I_u = setdiff(1:N_s, I_a);
I_a_bkw = sort(N_s - I_a + 1); 
%the indices of the available samples after time flip
I_u_bkw = sort(N_s - I_u + 1); 
%the indices of the unavailable samples after time flip
iter=0;
while 1
iter = iter+1;
U = zeros(length(I_u));
V = zeros(length(I_u),1);
for k = 1:N_omg
Hmtx = convmtx(H(:,k)', N_s-M+1);
Ak = Hmtx(:,I_a_bkw);
Uk = Hmtx(:, I_u_bkw);
U = U + Uk' * Uk;
mu_k = beta(k) * exp(ii*omega(k)*(N_s:-1:M)).';
V = V + Uk' * ( mu_k - Ak*flipud(sig(I_a)) );
end
y_u = U\V;
sig(I_u) = flipud(y_u);
%=================================================
%Either of the two lines below can be used.
[H, beta, omega, ~] = APES(sig, M, N_omg);
% [H, beta, omega, ~] = APES_FFT(sig, M, N_omg);
%=================================================
err_i1 = error_APES(sig, H, beta, omega);
if abs(err_i1 - err_i0)/abs(err_i0) < 1e-3
break;
else
disp(['iter: ', num2str(iter), ', error: ' num2str(err_i1)])
err_i0 = err_i1;
end
end
sig_out = sig;
end

