function [] = PlotDopplerSNR(SNR_db, v_mean_e, v_spread_e, OI, PI, Omega_rpm, mean_Phi)
    Phi_interest = mean_Phi(PI) .* 180/pi;
    Omega_interest = Omega_rpm(OI);
    figure; plot(SNR_db, squeeze(v_mean_e(:, OI, PI)), 'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2);
    grid on;
    xlabel('SNR [dB]', 'FontSize', 16)
    ylabel('Error in Mean Doppler velocity [m/s]', 'FontSize', 16);
    title(['Error in Mean Doppler velocity at ', num2str(Omega_interest), 'RPM', ', \Phi = ', num2str(Phi_interest), ' deg'], 'FontSize', 16);
    
    
    figure; plot(SNR_db, squeeze(v_spread_e(:, OI, PI)), 'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2);
    grid on;
    xlabel('SNR [dB]', 'FontSize', 16)
    ylabel('Error in Doppler spectrum width [m/s]', 'FontSize', 16);
    
    title(['Error in Doppler spectrum width at ', num2str(Omega_interest), 'RPM', ', \Phi = ', num2str(Phi_interest), ' deg'], 'FontSize', 16);

end