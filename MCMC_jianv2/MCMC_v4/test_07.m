%%
clear;

lambda = 3e-2;
dt = 1e-3;

mu = 5*4*pi/lambda*dt;
sig = 1*4*pi/lambda*dt;
Nu = 10000;

u = linspace(0, 15, Nu)*4*pi/lambda*dt+1e-7;


U = normrnd(mu, sig, [1 Nu])*4*pi/lambda*dt+1e-7;

[umn1, umn2] = meshgrid(U, U);

K = [0:1:4 101:1:105 201:1:205 301:1:305];

Nt = length(K);

Ns = 10000;

Np = 5;
M = 2;
Nscan = 100;

% 
% T1 = - 2 * Ns * M * Np;
% xmn = (umn1 - umn2);
% 
% T2i = -sin(Np .* xmn/2)./sin(xmn./2) .* ...
%     sin(Nscan * M .* xmn/2)./sin(Nscan .* xmn./2) .* cos(xmn/2 .* ((Np - 1) + (M - 1)*Nscan));
% T2 = sum(sum(T2i));


for i = 1:Nu
%     xm = (U - u(i));
%     T3i = 2 * sin(Np .* xm/2)./sin(xm./2) .* ...
%     sin(Nscan * M .* xm/2)./sin(Nscan .* xm./2) .* cos(xm/2 .* ((Np - 1) + (M - 1)*Nscan));
%     T3(i) = sum(T3i);
    F(i) = -sum(1 + Ns^2 * exp(-K.^2*sig^2) - 2 * Ns * exp(-K.^2.*sig^2/2) .* cos(K .* (mu - u(i))));
end

% F = T1 + T2 + T3 - max(T1+T2+T3);

figure(101); hold on; plot(u/(4*pi/lambda*dt), F-max(F)); grid on;



%% Approximation 


% A = Np .* M;
% % B = 1/24 * (Np * M * (2 * Np^2 + 3 * Np * (Nscan * (M - 1) - 1) + Nscan^2*(2*M^2 - 3*M+1) - 3*Nscan*(M-1)+1));
% % C = 1/5760 .* Np * M * (6 * Np^4 + 15 * Np^3 * (Nscan * (M - 1) - 1) + 10 * Np^2*(Nscan^2*(2*M^2 - 3*M + 1)...
% %     -3 * Nscan * (M - 1) + 1) + 15 * Np * Nscan * (M - 1) * (Nscan^2*(M - 1)*M - 2*Nscan*M + ...
% %     Nscan + 1) + Nscan^4*(6*M^4 - 15*M^3 + 10*M^2 - 1) - 15 * Nscan^3*(M-1)^2*M + ...
% %     5*Nscan^2*(2*M^2-3*M+1) - 1);
% 
% B = (1/2 * Nscan * (-(Np^3/(12 * Nscan)) + ...
%      1/2 * A * (Nscan/6 + (2 * (1/12 - 1/4 * (-1 + Np + Nscan * (-1 + M))^2))/Nscan)) * M - ...
%   1/24 * Np * Nscan^2 * M^3);
% 
% 
% C = (1/2 * Nscan * (Np^5/(960 * Nscan) - ...
%      1/48 * Np^3 * (Nscan/6 + (2 * (1/12 - 1/4 * (-1 + Np + Nscan * (-1 + M))^2))/Nscan) + ...
%      1/2 * Np * ((7 * Nscan^3)/1440 + ...
%         1/12 * Nscan * (1/12 - 1/4 * (-1 + Np + Nscan * (-1 + M))^2) + ...
%         (2 * (7/2880 - 1/96 * (-1 + Np + Nscan * (-1 + M))^2 + ...
%            1/192 * (-1 + Np + Nscan * (-1 + M))^4))/Nscan)) * M - ...
%   1/48 * Nscan^3 * (-(Np^3/(12 * Nscan)) + ...
%      1/2 * Np * (Nscan/6 + (2 * (1/12 - 1/4 * (-1 + Np + Nscan * (-1 + M))^2))/Nscan)) * M^3 + ...
%      (Np * Nscan^4 * M^5)/1920);
% 
% F2 = T1 - (A - B*sig^2+3*C*sig^4) - 2.*(A - B*((mu - u).^2 + sig.^2)+...
%     C/1e4.*((mu - u).^4 + 6 .* (mu - u).^2 .* sig.^2 + 3 .* sig.^4));
% 
% hold on; plot(u/(4*pi/lambda*dt), F2-max(F2)); grid on;
