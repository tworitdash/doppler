function [] = Plot2DDoppler(vel_axis, mean_Phi, Signal, SI, OI, SNR_db, Omega_rpm) 
    SNR_interest = SNR_db(SI);
    Omega_interest = Omega_rpm(OI);
    s = squeeze(mean(abs(squeeze(Signal(OI).doppler(:, SI, :, :))), 1));
    figure; imagesc(vel_axis, mean_Phi * 180/pi, abs(s)); shading flat; colorbar; %colormap('jet'); 
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Doppler spectrum', 'FontSize', 16);
    title(['Doppler spectrum at ', num2str(Omega_interest), 'RPM', ', SNR = ', num2str(SNR_interest), ' dB'], 'FontSize', 16);

end