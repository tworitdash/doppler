function [] = Plot2DDoppler(vel_axis, Phi, Signal, BI, SI, OI, BW_deg, SNR_db, Omega_rpm) 
    SNR_interest = SNR_db(SI);
    Omega_interest = Omega_rpm(OI);
    BW_interest = BW_deg(BI);
    s = squeeze(mean(abs(squeeze(Signal(OI, BI).doppler(:, SI, :, :))), 1));
    figure; imagesc(vel_axis, Phi(BI).Phi * 180/pi, abs(s)); shading flat; colorbar; %colormap('jet'); 
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Doppler spectrum', 'FontSize', 16);
    title(['Doppler spectrum at ', num2str(Omega_interest), 'RPM', ', SNR = ', num2str(SNR_interest), ' dB', ', BW = ', num2str(BW_interest), ' deg'], 'FontSize', 16);

end