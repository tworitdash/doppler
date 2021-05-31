clc
clear
close all

M0 = 1;
M1 = 5;
M2 = 1;
L = round(linspace(4,256,10));
vmax = 7.5;


for i = 1:numel(L)
    
    dv = 2*vmax/L(i);
    v1 = -vmax:dv:vmax;
    
    [X_PS_nf,M0_A(i),M1_A(i),M2_A(i)] = weather_sig_simulator_beta(M0,M1,M2,v1,dv);

    v2 = -vmax:dv:vmax-dv;
    [M0_B(i),M1_B(i),M2_B(i)] = gauss_calc_2(abs(X_PS_nf(1:end-1)),v2,dv);
    
end

M0_e = sqrt((M0_A-M0).^2);
M1_e = sqrt((M1_A-M1).^2);
M2_e = sqrt((M2_A-M2).^2);

M0_e_in = sqrt((M0_A-M0_B).^2);
M1_e_in = sqrt((M1_A-M1_B).^2);
M2_e_in = sqrt((M2_A-M2_B).^2);
    
figure
plot(L,(M0_e_in)),hold on
plot(L,(M1_e_in))
plot(L,(M2_e_in)),hold off
xlabel('Number of Doppler bins')
ylabel('Error')
legend('M0','M1','M2')
title('Error estimation for each moment between input and symmetric spectra')

figure
plot(L,(M0_e)),hold on
plot(L,(M1_e))
plot(L,(M2_e)),hold off
xlabel('Number of Doppler bins')
ylabel('Error')
legend('M0','M1','M2')
title('Error estimation for each moment when last bit is missing')

figure
plot(L,M0_A),hold on
plot(L,M0_B),hold off
title('M0 for both symmetric and un-symmetric spectra')
legend('Symmetric','Un-symmetric')

figure
plot(L,M1_A),hold on
plot(L,M1_B),hold off
title('M1 for both symmetric and un-symmetric spectra')
legend('Symmetric','Un-symmetric')

figure
plot(L,M2_A),hold on
plot(L,M2_B),hold off
title('M2 for both symmetric and un-symmetric spectra')
legend('Symmetric','Un-symmetric')

    
% figure
% plot(v1,(abs(X_PS_nf)))
% grid on
% grid minor
% title('Spectra with Simetric axis')

