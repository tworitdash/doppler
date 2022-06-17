clear;
close all;

%% Gibbs sampler for 2 variables

dt = 1;

omegamax = 1/(2*dt) * 2 * pi;

t = linspace(-256, 255, 512);

%% 
nt = length(t);

SNR_db = 0;
SNR = 10^(SNR_db/10);
omega1 = 0.3; omega2 = 0.31; b11 = 0.5403; b12 = -0.8415; b21 = -0.4161; b22 = -0.9093;
% omega1 = 0.3; omega2 = 0; b11 = 0.5403; b12 = -0.8415; b21 = 0; b22 = 0;
% s = cos(2 * pi * f * t);
s = b11 * cos(omega1 * t) + b12 * sin(omega1 * t) + b21 * cos(omega2 * t) + b22 * sin(omega2 * t);


sig_n2 = sum(abs(s).^2)/(nt * SNR);
z = s + sqrt(sig_n2) * (randn(1, nt));

%% MCMC with Gibbs sampling: 
%% Task: Estimate f and sig_n from the data with MCMC Gibbs sampling 

K = 10;
Ms = 1;
nu = nt - 1;

omega10 = 0.26;
omega20 = 0.26;

B110 = 0.7;
B120 = -0.6;
B210 = -0.5;
B220 = -1;

sig02 = sig_n2;


for m = 1:Ms
    disp(m);
    for k = 1:K
        D = z;
%% ========================================= For B11

        DB11 = z - B120 * sin(omega10 * t) - B210 * cos(omega20 * t) - B220 * sin(omega20 * t);
        B11khat = sum(DB11 .* cos(omega10.*t))/sum((cos(omega10.*t)).^2);
        X11 = cos(omega10 .* t).';
        cov_B11 = sig02 .* inv(X11.' * X11);
        
        B11(k, m) = normrnd(B11khat, sqrt(cov_B11), 1);
        
%% ========================================= For B12       
        DB12 = z - B11(k, m) * cos(omega10 * t) - B210 * cos(omega20 * t) - B220 * sin(omega20 * t);
        B12khat = sum(DB12 .* sin(omega10.*t))/sum((sin(omega10.*t)).^2);
        X12 = sin(omega10 .* t).';
        cov_B12 = sig02 .* inv(X12.' * X12);
        
        B12(k, m) = normrnd(B12khat, sqrt(cov_B12), 1); 
        
        
%% ========================================= For Omega1      
        
        omega01hat = omega10;
        
        for l = 2:50
            J = (-B11(k, m) .* t .* sin(omega01hat(l-1) .* t) + B12(k, m) .* t .* cos(omega01hat(l-1) .* t)).';
            func = (B11(k, m) .* cos(omega01hat(l-1) .* t) + B12(k, m) .* sin(omega01hat(l-1) .* t) +...
                B210 .* cos(omega20 .* t) + B220 .* sin(omega20 .* t)).';
            delD1 = D.' - func;
            delomega1 = inv(J.' * J) * J.' * delD1;
            omega01hat(l) = omega01hat(l-1)  + delomega1;
        end
        
        omega10 = omega01hat(end);
        
%         figure(100); hold on; plot(omega0hat/(2*pi));
        
        
        X_omega1 = (-t .* B11(k, m) .* sin(omega10*t) + t .* B12(k, m) .* cos(omega10*t)).';
%         Dhat = B1(k) .* cos(omega0.*t_a);
%         s_omega2 = 1/(na - 1) .* (D - Dhat).' * (D - Dhat);
        cov_Omega1 = sig02 * inv(X_omega1.' * X_omega1);
        
        Omega1(k, m) = normrnd(omega10, sqrt(cov_Omega1), 1);
        
%% ========================================= For B21    
        
        DB21 = z - B11(k, m) * cos(Omega1(k, m) * t) - B12(k, m) * sin(Omega1(k, m) * t) - B220 * sin(omega20 * t);
        B21khat = sum(DB21 .* cos(omega20.*t))/sum((cos(omega20.*t)).^2);
        X21 = cos(omega20 .* t).';
        cov_B21 = sig02 .* inv(X21.' * X21);
        
        B21(k, m) = normrnd(B21khat, sqrt(cov_B21), 1);
%% ========================================= For B22      

        DB22 = z - B11(k, m) * cos(Omega1(k, m) * t) - B12(k, m) * sin(Omega1(k, m) * t) - B21(k, m) * cos(omega20 * t);
        B22khat = sum(DB22 .* sin(omega20.*t))/sum((sin(omega20.*t)).^2);
        X22 = sin(omega20 .* t).';
        cov_B22 = sig02 .* inv(X22.' * X22);
        
        B22(k, m) = normrnd(B22khat, sqrt(cov_B22), 1);
        

        
%% ========================================= For Omega2
        omega02hat = omega20;
        
        for l = 2:50
            J = (-B21(k, m) .* t .* sin(omega02hat(l-1) .* t) + B22(k, m) .* t .* cos(omega02hat(l-1) .* t)).';
            func = (B11(k, m) .* cos(Omega1(k, m) .* t) + B12(k, m) .* sin(Omega1(k, m) .* t) +...
                B21(k, m) .* cos(omega02hat(l-1) .* t) + B22(k, m) .* sin(omega02hat(l-1) .* t)).';
            delD2 = D.' - func;
            delomega2 = inv(J.' * J) * J.' * delD2;
            omega02hat(l) = omega02hat(l-1)  + delomega2;
        end
        
        omega20 = omega02hat(end);
        
%         figure(100); hold on; plot(omega0hat/(2*pi));
        
        
        X_omega2 = (-t .* B21(k, m) .* sin(omega20*t) + t .* B22(k, m) .* cos(omega20*t)).';
%         Dhat = B1(k) .* cos(omega0.*t_a);
%         s_omega2 = 1/(na - 1) .* (D - Dhat).' * (D - Dhat);
        cov_Omega2 = sig02 * inv(X_omega2.' * X_omega2);
        
        Omega2(k, m) = normrnd(omega20, sqrt(cov_Omega2), 1);
        
        
%% Initializing for the next iteration 
        
        
        omega10 = Omega1(k, m);
        omega20 = Omega2(k, m);
        B110 = B11(k, m);
        B120 = B12(k, m);
        B210 = B21(k, m);
        B220 = B22(k, m);
%         freq(k, m) = mod(Omega(k, m)/(2*pi), fmax);
%         vel(k, m) = mod(freq(k, m)*lambda/2, fmax*lambda/2);
    end
end
% figure; plot(B1(:, 1)); figure; plot(vel(:, 1));
figure; histogram(B11(:, 1)); title('b_{11}');
hold on; histogram(b11, 100); 

figure; histogram(B12(:, 1)); title('b_{12}');
hold on; histogram(b12, 100); 

figure; histogram(B21(:, 1)); title('b_{21}');
hold on; histogram(b21, 100); 

figure; histogram(B22(:, 1)); title('b_{22}');
hold on; histogram(b22, 100); 

figure; histogram(Omega1(:, 1)); title('\omega_1');
hold on; histogram(omega1, 100);

figure; histogram(Omega2(:, 1));title('\omega_2');
hold on; histogram(omega2, 100);