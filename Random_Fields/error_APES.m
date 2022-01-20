function error = error_APES(sig, H, beta, omega)
ii = sqrt(-1);
if size(sig,1)==1
sig = sig.';
end
N_s = length(sig);
[m, N_omg] = size(H);
error = 0;
for k = 1: N_omg
Hmtx = convmtx(H(:,k)', N_s-m+1);
error0 = Hmtx * sig(end:-1:1) - beta(k) * exp(ii*omega(k)*(N_s:-1:m)).';
error = error + norm(error0)^2;
end
end
