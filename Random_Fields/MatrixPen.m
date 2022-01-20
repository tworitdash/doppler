%% Signal model with HD
clear;
close all;

markers = load('../mono/markers.mat');
markers = markers.markers;

N = 1;
r0 = 100;
dr = 50;
ph0 = 0;
dph = 1.8*pi/180;
n_rot = 2;
N_Sweep = 11;
N_sec = 100;
Nt = n_rot*N_sec*N_Sweep;
SNR_db = 30;
lambda = 0.03;
dt = 1e-3;
u_mean = 1;
u_sigma = 0.1;
v_mean = 0;
v_sigma = 0;

[z_model] = Zmodel(N, r0, dr, ph0, dph, Nt, SNR_db, lambda, dt, u_mean, u_sigma, v_mean, v_sigma); 

[ZFFT, PT, mu, sigma, vel_axis, dv] = Spec(z_model, Nt, dt, lambda, SNR_db, 1, 1, 2);

%% Decide the samples available and club them

[zi, z] = Zavail(z_model, n_rot, N_Sweep, N_sec);

[~,~, mu_avail, sigma_avail, vel_axis_avail, dv_avail] = Spec(zi(1, :), N_Sweep, dt, lambda, SNR_db, 1, 1, 6);
%% Finding signal order of the segments of the signal defined with n_rot 

L = 5;
M = N_Sweep - L + 1;
for i = 1:n_rot
    for l = 1:L
        H(i, l, :) = zi(i, l:l+M-1);
    end
    [U(i, :, :), S(i, :, :), V(i, :, :)] = svd(squeeze(H(i, :, :)));
    Uktb(i, :, :) = [squeeze(U(i, 2:end, :)) squeeze(U(i, 1:end-1, :))];
end

figure; plot((abs(diag(squeeze(S(1, :, :))))));

%% MP Signal Reconstruction STEP 1

Mhat = 4;

M1 = size(zi, 2); M2 = size(zi, 2);
Mi = [M1 M2];

M1M2min = min(M1 - Mhat, M2 - Mhat);
Lmp = round(1./2 .* (Mhat + M1M2min));

for i = 1:n_rot
    for k = 1:Lmp+1
        D(i, k, :) = zi(i, k:Mi(i)-Lmp+k-1);
    end
end

for i = 1:n_rot
    for k = 1:2
        Han(i, k, :, :) = [squeeze(D(i, k:Lmp+k-1, :))];
    end
end

for i = 1:n_rot
    X(i, :, :) = [squeeze(Han(1, 1, :, :)); squeeze(Han(2, 1, :, :))];
    [Ump(i, :, :), Smp(i, :, :), Vmp(i, :, :)] = svd(squeeze(X(i, :, :)));
    Xt(i, :, :) = squeeze(Ump(i, 1:Mhat, 1:Mhat)) * squeeze(Smp(i, 1:Mhat, 1:Mhat)) * squeeze(Vmp(i, 1:Mhat, 1:Mhat))';
end

S1 = squeeze(Smp(1, 1:Mhat, 1:Mhat));
S2 = squeeze(Smp(2, 1:Mhat, 1:Mhat));

U1 = squeeze(Ump(1, 1:Mhat, 1:Mhat));
U2 = squeeze(Ump(2, 1:Mhat, 1:Mhat));

V1 = squeeze(Vmp(1, 1:Mhat, 1:Mhat));
V2 = squeeze(Vmp(2, 1:Mhat, 1:Mhat));

A = inv(S1) * U1' * U2 * S2 * V2' * V1;

ei = eig(A); 

s = [zi(1, :) zi(2, :)];

for i = 1:n_rot
    Nvec((i - 1)*N_Sweep+1:i*N_Sweep) = N_sec*N_Sweep*(i - 1)+1:N_sec*N_Sweep*(i - 1)+N_Sweep;
end

for l = 1:length(Nvec)
    Z(l, :) = ei.^(Nvec(l) - 1);
end

acoeff = pinv(Z) * s.';

%% Signal Reconstruction STEP 2
valid = 1;
count = 1; epi(1) = 0;

while valid

count = count + 1;

for m = 1:Nt
    shat(m) = sum(acoeff .* ei.^(m-1));
end

for i = 1:n_rot
    sihat(i, :) = shat(Nvec((i - 1)*N_Sweep+1:i*N_Sweep));
end

% Error second norm

epi(count) = sqrt(sum((sihat(1, :) - zi(1, :)).^2)) + sqrt(sum((sihat(2, :) - zi(2, :)).^2));


if (epi(count) <= epi(count-1)) && (count > 2)
    valid = 0;
end

%%  Signal Reconstruction STEP 3

shat(Nvec) = [zi(1, :) zi(2, :)];

Lmp1 = round((Nt/3+Nt/2)/2);

for l = 1:Nt-Lmp1
    Y(l, :) = shat(l:Lmp1+l);
end

[Uy, Sy, Vy] = svd(Y);

% figure; plot((abs(diag(Sy))));

S_norm = abs(diag(Sy)).^2./max(abs(diag(Sy)).^2);

[~, idxs] = find(S_norm > 1e-3);

Trunc = length(idxs);

Uyt = Uy;
Vyt = Vy(1:Trunc, 1:Trunc);
Syt = Sy(:, 1:Trunc);

Y1 = Uyt * Syt * Vyt(1:end-1, :)';
Y2 = Uyt * Syt * Vyt(2:end, :)';

Ay = pinv(Y1)*Y2;

clear ei; 
clear acoeff;


ei = eig(Ay);


for l = 1:Nt
    Zy(l, :) = ei.^(l - 1);
end

acoeff = pinv(Zy) * shat.';
clear Zy;

end

%% Comparison with zero padding

Z_model_zero = z_model;
Z_model_zero(Nvec) = z_model(Nvec);
Z_model_zero(setdiff(1:end, Nvec)) = 0;

Spec(shat, length(shat), dt, lambda, SNR_db, 1, 1, 10);
Spec(Z_model_zero, length(Z_model_zero), dt, lambda, SNR_db, 1, 3, 12);


% figure(108); hold on; plot(vel_axis, db(abs(fftshift(fft(shat)))));
% figure(108); hold on; plot(vel_axis, db(abs(fftshift(fft(Z_model_zero)))));