close all;
clear;
clc;

f = 10e9; % place the frequencies used in FEKO.
c = 3e8;
lamb = c./f;

% close all;


% T = readtable('/tudelft.net/staff-groups/ewi/me/MS3/Studenten/Tworit Dash/Extra Project/MicroDoppler/2Blades_3D/Data/E_field_theta_amp_HV_2Blades.dat');
% 
% T1 = readtable('/tudelft.net/staff-groups/ewi/MS3/Studenten/Tworit Dash/Extra Project/MicroDoppler/2Blades_3D/Data/E_field_theta_phase_HV_2Blades.dat');

% T2 = readtable('/tudelft.net/staff-groups/ewi/me/MS3/Studenten/Tworit Dash/Extra Project/MicroDoppler/1Blade_3D/Data/E_field_phi_amp_HH.dat');
% T3 = readtable('/tudelft.net/staff-groups/ewi/me/MS3/Studenten/Tworit Dash/Extra Project/MicroDoppler/1Blade_3D/Data/E_field_phi_phase_HH.dat');

% All er = 5

T = readtable('Cylinder_Doppler_SP.dat');


% % eth_amp = T.FarField1_dBV_;
% eph_amp = T2.x1Bladde3D_all_er5_dBV_;
% % eth_phase = T1.FarField1_deg_;
% eph_phase = T3.x1Bladde3D_all_er5_deg_;

eph_amp = T.Magnitude_dBV_;
eph_phase = T.Phase_deg_;


% eth = 10.^(eth_amp./20) .* exp(1j .* eth_phase .* pi/180);
eph = 10.^(eph_amp./20) .* exp(1j .* eph_phase .* pi/180);


%Micro Doppler

wlen = 32;                        % window length (recomended to be power of 2)
hop = 31;
% hop = 64; % hop size (recomended to be power of 2)
nfft = wlen;                        % number of fft points (recomended to be power of 2)

% Calculating Fs from PRI
Omega = 60000; % rpm
Omega_s = Omega/60; % rps

% L = 21.7e-2+2.5e-2; %length of the blades
L = 38;
V_tip = 2 * pi * L * Omega_s;

f_dmax = (2 * V_tip)./lamb;

fs =  Omega_s * 360;
% fs = 1e6;


[S1,F1,Ti1,P1] = spectrogram(eph, wlen, hop, nfft, fs); % For HH


F1_doppler = F1 - fs/2; % Doppler frequency [Hz]


S1 = S1./max(S1(:));
S1_norm_db = 20*log(abs(fftshift(S1',2)));

% i = i + 1;
figure;
imagesc(F1_doppler*lamb/2, Ti1, S1_norm_db);
colormap('jet');
colorbar;

xlabel('Dopplervelocity [ms^{-1}]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Time[s]', 'FontSize', 12, 'FontWeight', 'bold');
title([' Micro Doppler N_{fft} = ', num2str(nfft), ', overlap = ', num2str(hop),...
    ', fs = ', num2str(fs*1e-6), ' MHz, ', 'WINDOW = ', num2str(wlen), ' Material = PEC'], 'FontSize', 10, 'FontWeight', 'bold');
% title([' Micro Doppler N_{fft} = ', num2str(nfft), ', overlap = ', num2str(hop),...
%     ', fs = ', num2str(fs*1e-6), ' MHz, ', 'WINDOW = ', num2str(wlen), ', For side ', s], 'FontSize', 10, 'FontWeight', 'bold');

%xlim([-50 50]);



